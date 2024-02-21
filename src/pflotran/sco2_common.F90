module SCO2_Common_module

#include "petsc/finclude/petscsys.h"

    use SCO2_Aux_module
    use Global_Aux_module
    use PFLOTRAN_Constants_module
    use petscsys

    implicit none

    PetscReal, parameter :: eps = 1.d-8
    PetscReal, parameter :: floweps = 1.d-24

    public :: SCO2Accumulation, &
              SCO2Flux, &
              SCO2BCFlux, &
              SCO2AuxVarComputeAndSrcSink, &
              SCO2AccumDerivative, &
              SCO2FluxDerivative, &
              SCO2BCFluxDerivative, &
              SCO2SrcSinkDerivative

contains

! ************************************************************************** !

subroutine SCO2Accumulation(sco2_auxvar,global_auxvar,material_auxvar, &
                            soil_heat_capacity,option,Res)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof)
  PetscInt :: iphase, icomp
  PetscReal :: porosity
  PetscReal :: volume_over_dt

  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  porosity = sco2_auxvar%effective_porosity

  ! accumulation term units = kg/s
  Res = 0.d0
  ! Trapped gas should be accounted for in the gas phase
  do iphase = 1, option%nphase - 1
    ! Res[kg comp/sec] =      sat[m^3 phase/m^3 void] *
    !                         den_kg[kg phase/m^3 phase] *
    !                         xmass[kg comp/kg phase] *
    !                         por[m^3 void/m^3 bulk] *
    !                         vol/dt[m^3 bulk/sec]
    do icomp = 1, option%nflowspec - 1
      Res(icomp) = Res(icomp) + ( sco2_auxvar%sat(iphase) * &
                            sco2_auxvar%den_kg(iphase) * &
                            sco2_auxvar%xmass(icomp,iphase) ) * &
                            porosity * volume_over_dt
    enddo
  enddo

  ! Salt precipitate is calculated as fraction of total porosity,
  ! so it needs to be added separately.
  Res(SCO2_SALT_EQUATION_INDEX) = Res(SCO2_SALT_EQUATION_INDEX) + &
                              sco2_auxvar%m_salt(TWO_INTEGER) * &
                              volume_over_dt

  if (.not. sco2_isothermal) then
    do iphase = 1, option%nphase
      ! Res[MJ/s] =    sat[m^3 phase/m^3 void] *
      !                    den_kg[kg phase/m^3 phase] * U[MJ/kg phase] *
      !                    por[m^3 void/m^3 bulk] *
      !                    vol/dt[m^3 bulk/sec]
      Res(SCO2_ENERGY_EQUATION_INDEX) = Res(SCO2_ENERGY_EQUATION_INDEX) + &
                                        sco2_auxvar%sat(iphase) * &
                                        sco2_auxvar%den_kg(iphase) * &
                                        sco2_auxvar%U(iphase) * &
                                        porosity * volume_over_dt
    enddo
    ! Add rock component
    Res(SCO2_ENERGY_EQUATION_INDEX) = Res(SCO2_ENERGY_EQUATION_INDEX) + &
                                      (1.d0 - porosity) * &
                                      material_auxvar%soil_particle_density * &
                                      soil_heat_capacity * sco2_auxvar%temp * &
                                      volume_over_dt
  endif
end subroutine SCO2Accumulation

! ************************************************************************** !

subroutine SCO2Flux(sco2_auxvar_up,global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_cc_up, &
                    sco2_auxvar_dn,global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_cc_dn, &
                    area, dist, upwind_direction_, &
                    option,v_darcy,Res,&
                    update_upwind_direction_, &
                    count_upwind_direction_flip_)
  !
  ! Computes the internal flux terms for the residual
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use Material_Aux_module
  use Connection_module
  use Fracture_module
  use Upwind_Direction_module
  use Characteristic_Curves_Thermal_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar_up, sco2_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  class(cc_thermal_type) :: thermal_cc_up, thermal_cc_dn
  PetscReal :: Res(option%nflowdof)
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: iphase, icomp

  PetscInt :: lid, gid, pid, tgid, wid, co2_id, sid
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: delta_temp
  PetscReal :: uH
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: gravity_term
  PetscReal :: mobility, q
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: dkeff_up_dsatlup, dkeff_up_dTup, dkeff_dn_dsatldn, dkeff_dn_dTdn

  PetscReal :: xmass(option%nflowspec)
  PetscReal :: tot_mass_flux, component_mass_flux, co2_mass_flux, &
               co2_mole_flux, salt_mass_flux, &
               water_mass_flux, salt_diff_flux
  PetscReal :: delta_xmass, delta_xmol, den_dn, den_up, density_ave
  PetscReal :: den_kg_up, den_kg_dn, density_kg_ave
  PetscReal :: sat_dn, sat_up
  PetscReal :: stpd_ave_over_dist, stpd_up, stpd_dn
  PetscReal :: al, alp
  PetscReal :: dheat_flux_ddelta_temp
  PetscReal :: dtot_mole_flux_ddeltaX
  PetscReal :: dsalt_mass_flux_ddeltaX
  PetscReal :: up_scale, dn_scale
  PetscBool :: upwind
  PetscReal :: visc_mean, kr
  PetscReal :: tempreal
  PetscReal, parameter :: epsilon = 1.d-20

  lid = option%liquid_phase
  gid = option%gas_phase
  pid = option%precipitate_phase
  tgid = option%trapped_gas_phase

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  perm_up = perm_up * sco2_auxvar_up%effective_permeability
  perm_dn = perm_dn * sco2_auxvar_dn%effective_permeability

  ! Harmonic permeability
  perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)

  Res = 0.d0

  v_darcy = 0.d0

  do iphase = 1 , option%nphase - 1 ! No advection or diffusion through salt phase

    if (sco2_auxvar_up%sat(iphase) == 0.d0 .and. &
        sco2_auxvar_dn%sat(iphase) == 0.d0) cycle

    ! Harmonic mean on viscosity
    visc_mean = (sco2_auxvar_up%visc(iphase) * sco2_auxvar_dn%visc(iphase) * &
               (dist_up + dist_dn)) / (sco2_auxvar_up%visc(iphase) * dist_up + &
                sco2_auxvar_dn%visc(iphase) * dist_dn)
    ! STOMP takes harmonic mean on density when old velocity is 0

    ! Advection

    if (sco2_auxvar_up%mobility(iphase) + &
        sco2_auxvar_dn%mobility(iphase) > eps) then

      density_kg_ave = SCO2AverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          sco2_auxvar_up%den_kg, &
                                          sco2_auxvar_dn%den_kg)

      gravity_term = density_kg_ave * dist_gravity
      delta_pressure = sco2_auxvar_up%pres(iphase) - &
                       sco2_auxvar_dn%pres(iphase) + &
                       gravity_term

      up_scale = 0.d0
      dn_scale = 0.d0
      upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                               PETSC_FALSE, &
                               count_upwind_direction_flip_, &
                               liq_upwind_flip_count_by_res, &
                               liq_upwind_flip_count_by_jac)
      if (upwind) then
        up_scale = 1.d0
        mobility = sco2_auxvar_up%mobility(iphase)
        kr = sco2_auxvar_up%kr(iphase)
        xmass(:) = sco2_auxvar_up%xmass(:,iphase)
        uH = sco2_auxvar_up%H(iphase)
        density_kg_ave = sco2_auxvar_up%den_kg(iphase)
        !perm_ave_over_dist = perm_up / (dist_up + dist_dn)
      else
        dn_scale = 1.d0
        mobility = sco2_auxvar_dn%mobility(iphase)
        kr = sco2_auxvar_dn%kr(iphase)
        xmass(:) = sco2_auxvar_dn%xmass(:,iphase)
        uH = sco2_auxvar_dn%H(iphase)
        density_kg_ave = sco2_auxvar_dn%den_kg(iphase)
        !perm_ave_over_dist = perm_dn / (dist_up + dist_dn)
      endif

      mobility = kr / visc_mean

      if (mobility > floweps .and. dabs(delta_pressure) > 0.d0 ) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * &
                          delta_pressure
        ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area
        ! mass_flux[kg phase/sec] = q[m^3 phase/sec] *
        !                             density_ave[kg phase/m^3 phase]
        tot_mass_flux = q * density_kg_ave
        ! comp_mass_flux[kg comp/sec] = tot_mass_flux[kg phase/sec] *
        !                                 xmass[kg comp/kg phase]

        do icomp = 1 , option%nflowspec - 1 !Handle salt separately

          component_mass_flux = tot_mass_flux * xmass(icomp)
          Res(icomp) = Res(icomp) + component_mass_flux

        enddo
        if (.not. sco2_isothermal) then
          ! Energy flux
          Res(SCO2_ENERGY_EQUATION_INDEX) = Res(SCO2_ENERGY_EQUATION_INDEX) + &
                                            tot_mass_flux * uH
        endif
      endif
    endif

    ! Diffusion

    ! Compute mole flux for CO2, mass flux for NaCl

    ! Harmonic diffusion coefficient

    if (sco2_harmonic_diff_density) then
        density_ave = 1.d0
        den_up = sco2_auxvar_up%den(iphase)
        den_dn = sco2_auxvar_dn%den(iphase)
        den_kg_up = sco2_auxvar_up%den_kg(iphase)
        den_kg_dn = sco2_auxvar_dn%den_kg(iphase)
    else
      den_up = 1.d0
      den_dn = 1.d0
      den_kg_up = 1.d0
      den_kg_dn = 1.d0
      ! use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = SCO2AverageDensity(iphase, &
                                        global_auxvar_up%istate, &
                                        global_auxvar_dn%istate, &
                                        sco2_auxvar_up%den, &
                                        sco2_auxvar_dn%den)
      density_kg_ave = SCO2AverageDensity(iphase, &
                                        global_auxvar_up%istate, &
                                        global_auxvar_dn%istate, &
                                        sco2_auxvar_up%den_kg, &
                                        sco2_auxvar_dn%den_kg)
    endif

    if (iphase == LIQUID_PHASE) then
      ! CO2 Mole Flux in the aqueous phase
      ! Include diffusion and longitudinal dispersion
      stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(co2_id,iphase) + &
                 sco2_auxvar_up%dispersivity(co2_id,iphase) * &
                 v_darcy(iphase))
      stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(co2_id,iphase) + &
                 sco2_auxvar_dn%dispersivity(co2_id,iphase) * &
                 v_darcy(iphase))

      ! units = [kg/m^2/s bulk]
      stpd_ave_over_dist = stpd_up*stpd_dn / &
                           (stpd_up*dist_dn + stpd_dn*dist_up)

      ! units = kg/sec
      dtot_mole_flux_ddeltaX = stpd_ave_over_dist * area

      delta_xmol = sco2_auxvar_up%xmol(co2_id,iphase) * den_up - &
                   sco2_auxvar_dn%xmol(co2_id,iphase) * den_dn

      co2_mole_flux = dtot_mole_flux_ddeltaX * delta_xmol
    else
      ! Vapor Mole Flux in the gas phase
      ! Include diffusion and longitudinal dispersion
      stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(wid,iphase) + &
                 sco2_auxvar_up%dispersivity(wid,iphase) * &
                 v_darcy(iphase))
      stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(wid,iphase) + &
                 sco2_auxvar_dn%dispersivity(wid,iphase) * &
                 v_darcy(iphase))

      ! units = [kg/m^2/s bulk]
      stpd_ave_over_dist = stpd_up*stpd_dn / &
                           (stpd_up*dist_dn + stpd_dn*dist_up)

      ! units = kg/sec
      dtot_mole_flux_ddeltaX = stpd_ave_over_dist * area

      delta_xmol = sco2_auxvar_up%xmol(wid,iphase) * den_up - &
                   sco2_auxvar_dn%xmol(wid,iphase) * den_dn

      co2_mole_flux = -dtot_mole_flux_ddeltaX * delta_xmol
    endif

    ! Salt mass flux
    ! Patankar salt transport
    ! Include diffusion and longitudinal dispersion
    stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(sid,iphase) + &
               sco2_auxvar_up%dispersivity(sid,iphase) * &
               v_darcy(iphase))
    stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(sid,iphase) + &
               sco2_auxvar_dn%dispersivity(sid,iphase) * &
               v_darcy(iphase))

    tempreal = stpd_up*dist_up+stpd_dn*dist_dn
    if (tempreal > 0.d0) then
      stpd_ave_over_dist = stpd_up*stpd_dn / &
                           (stpd_up*dist_dn + stpd_dn*dist_up)
    else
      stpd_ave_over_dist = 0.d0
    endif

    al = max(v_darcy(iphase),0.d0) + stpd_ave_over_dist * max((1.d0 - &
         (1.d-1 * dabs(v_darcy(iphase))/(stpd_ave_over_dist + epsilon))) ** 5, &
          0.d0)
    alp = max(-v_darcy(iphase),0.d0) + stpd_ave_over_dist * max((1.d0 - &
         (1.d-1 * dabs(v_darcy(iphase))/(stpd_ave_over_dist + epsilon))) ** 5, &
          0.d0)

    salt_mass_flux = (al * sco2_auxvar_up%xmass(sid,iphase) * &
                          sco2_auxvar_up%den_kg(iphase) - &
                          alp * sco2_auxvar_dn%xmass(sid,iphase) * &
                          sco2_auxvar_dn%den_kg(iphase)) * area

    ! Diffusive component of salt flux
    ! units = kg/sec
    dsalt_mass_flux_ddeltaX = stpd_ave_over_dist * area

    delta_xmass =  sco2_auxvar_up%xmass(sid,iphase) * &
                   sco2_auxvar_up%den_kg(iphase) - &
                   sco2_auxvar_dn%xmass(sid,iphase) * &
                   sco2_auxvar_dn%den_kg(iphase)

    salt_diff_flux = dsalt_mass_flux_ddeltaX * delta_xmass

    if (iphase == ONE_INTEGER) then
        water_mass_flux = -1.d0 * fmw_comp(1) * &
                      (co2_mole_flux  + salt_diff_flux / fmw_comp(3))
    else
        water_mass_flux = -1.d0 * fmw_comp(1) * co2_mole_flux
    endif

    co2_mass_flux = co2_mole_flux * fmw_comp(2)

    ! Diffusive Contributions
    Res(SCO2_WATER_EQUATION_INDEX) = Res(SCO2_WATER_EQUATION_INDEX) + &
                                     water_mass_flux
    Res(SCO2_CO2_EQUATION_INDEX) = Res(SCO2_CO2_EQUATION_INDEX) + &
                                   co2_mass_flux
    Res(SCO2_SALT_EQUATION_INDEX) = Res(SCO2_SALT_EQUATION_INDEX) + &
                                    salt_mass_flux

    ! ! MAN: an effective multiphase diffusion coefficient approach:
    ! ! For CO2:
    ! do iphase = 1 , option%nphase - 1
    !   stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(co2_id,iphase) + &
    !              sco2_auxvar_up%dispersivity(co2_id,iphase) * &
    !              v_darcy(iphase))
    !   stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(co2_id,iphase) + &
    !              sco2_auxvar_dn%dispersivity(co2_id,iphase) * &
    !              v_darcy(iphase))

    !   ! Take the harmonic mean / dist:
    !   stpd_ave_over_dist = stpd_up*stpd_dn / &
    !                        (stpd_up*dist_dn + stpd_dn*dist_up)

    !   sigma(iphase) = stpd_ave_over_dist * area
    ! enddo

    ! multiphase_grad =(sco2_auxvar_up%xmol(co2_id,option%gas_phase) * den_up - &
    !               sco2_auxvar_dn%xmol(co2_id,option%gas_phase) * den_dn) / &
    !               (sco2_auxvar_up%xmol(co2_id,option%liquid_phase) * den_up - &
    !               sco2_auxvar_dn%xmol(co2_id,option%liquid_phase) * den_dn)
    ! co2_mole_flux = (sigma(ONE_INTEGER) + &
    !               sigma(TWO_INTEGER) * multiphase_grad) * &
    !               (sco2_auxvar_up%xmol(co2_id,option%liquid_phase) * den_up - &
    !               sco2_auxvar_dn%xmol(co2_id,option%liquid_phase) * den_dn)

    ! ! For salt:
    ! do iphase = 1 , option%nphase - 1
    !   stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(sid,iphase) + &
    !              sco2_auxvar_up%dispersivity(sid,iphase) * &
    !              v_darcy(iphase))
    !   stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(sid,iphase) + &
    !              sco2_auxvar_dn%dispersivity(sid,iphase) * &
    !              v_darcy(iphase))

    !   ! Take the harmonic mean / dist:
    !   stpd_ave_over_dist = stpd_up*stpd_dn / &
    !                        (stpd_up*dist_dn + stpd_dn*dist_up)

    !   sigma(iphase) = stpd_ave_over_dist * area
    ! enddo

    ! multiphase_grad =(sco2_auxvar_up%xmol(sid,option%gas_phase) * den_up - &
    !               sco2_auxvar_dn%xmol(sid,option%gas_phase) * den_dn) / &
    !               (sco2_auxvar_up%xmol(sid,option%liquid_phase) * den_up - &
    !               sco2_auxvar_dn%xmol(sid,option%liquid_phase) * den_dn)
    ! salt_mole_flux = (sigma(ONE_INTEGER) + &
    !               sigma(TWO_INTEGER) * multiphase_grad) * &
    !               (sco2_auxvar_up%xmol(sid,option%liquid_phase) * den_up - &
    !               sco2_auxvar_dn%xmol(sid,option%liquid_phase) * den_dn)


    ! water_mass_flux = -1.d0 * fmw_comp(1) * &
    !                   (co2_mole_flux + salt_mole_flux)

  enddo

  ! Conduction
  ! MAN: Need to extend the thermal conductivity functionality to include
  !      salt
  if (.not. sco2_isothermal) then
    sat_up = sco2_auxvar_up%sat(lid)
    sat_dn = sco2_auxvar_dn%sat(lid)

    ! derive wet and dry conductivities with anisotropy tensor and direction
    call thermal_cc_up%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)

    call thermal_cc_dn%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)

    ! thermal conductivity a function of temperature and liquid saturation
    call thermal_cc_up%thermal_conductivity_function%CalculateTCond(sat_up, &
         sco2_auxvar_up%temp,sco2_auxvar_up%effective_porosity, &
         k_eff_up,dkeff_up_dsatlup,dkeff_up_dTup,option)

    call thermal_cc_dn%thermal_conductivity_function%CalculateTCond(sat_dn, &
         sco2_auxvar_dn%temp,sco2_auxvar_dn%effective_porosity, &
         k_eff_dn,dkeff_dn_dsatldn,dkeff_dn_dTdn,option)

    if (k_eff_up > 0.d0 .or. k_eff_dn > 0.d0) then
      tempreal = k_eff_up*dist_dn + k_eff_dn*dist_up
      k_eff_ave = k_eff_up*k_eff_dn/tempreal
    else
      k_eff_ave = 0.d0
    endif

    ! units:
    ! k_eff = W/K-m = J/s/K-m
    ! delta_temp = K
    ! area = m^2
    ! heat_flux = k_eff * delta_temp * area = J/s
    ! 1.0E-6 term accounts for change in units: J/s -> MJ/s

    delta_temp = sco2_auxvar_up%temp - sco2_auxvar_dn%temp
    dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
    heat_flux = dheat_flux_ddelta_temp * delta_temp

    ! MJ/s or MW
    Res(SCO2_ENERGY_EQUATION_INDEX) = Res(SCO2_ENERGY_EQUATION_INDEX) + heat_flux
  endif

end subroutine SCO2Flux

! ************************************************************************** !

subroutine SCO2BCFlux(ibndtype, auxvar_mapping, auxvars, sco2_auxvar_up, &
                      global_auxvar_up, sco2_auxvar_dn, global_auxvar_dn, &
                      material_auxvar_dn, thermal_cc_dn, area, dist, &
                      upwind_direction_, option, v_darcy, &
                      Res, update_upwind_direction_, &
                      count_upwind_direction_flip_)
  !
  ! Computes boundary flux terms for the residual.
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use Material_Aux_module
  use Upwind_Direction_module
  use Characteristic_Curves_Thermal_module
  use Utility_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(16)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(sco2_auxvar_type) :: sco2_auxvar_up, sco2_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: Res(1:option%nflowdof)
  class(cc_thermal_type) :: thermal_cc_dn
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscInt :: iphase, icomp

  PetscInt :: lid, gid, pid, tgid, wid, co2_id, sid
  PetscReal :: perm_dn
  PetscReal :: delta_pressure
  PetscReal :: delta_temp
  PetscReal :: uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: gravity_term
  PetscReal :: mobility, q
  PetscReal :: kr, visc_mean
  PetscReal :: k_eff_dn, k_eff_ave, heat_flux

  PetscInt :: bc_type
  PetscReal :: boundary_pressure
  PetscReal :: xmass(option%nflowspec)
  PetscReal :: tot_mass_flux, component_mass_flux
  PetscReal :: co2_mass_flux, co2_mole_flux, salt_mass_flux, &
               salt_diff_flux
  PetscReal :: sat_dn
  PetscReal :: dn_scale
  PetscBool :: upwind
  PetscReal :: delta_xmol, den_dn, den_up
  PetscReal :: delta_xmass, den_kg_dn, den_kg_up, density_kg_ave
  PetscReal :: dheat_flux_ddelta_temp, dkeff_dn_dsatldn, &
               dkeff_dn_dtdn
  PetscReal :: al, alp
  PetscReal :: dsalt_mass_flux_ddeltax, &
               dtot_mole_flux_ddeltax, dv_darcy_ddelta_pressure
  PetscReal :: stpd_ave_over_dist, stpd_dn, stpd_up
  PetscReal :: dist_up, dist_dn
  PetscReal :: water_mass_flux
  PetscReal :: tempreal
  PetscInt :: idof
  PetscReal, parameter :: epsilon = 1.d-20

  lid = option%liquid_phase
  gid = option%gas_phase
  pid = option%precipitate_phase
  tgid = option%trapped_gas_phase

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  perm_dn = perm_dn * sco2_auxvar_dn%effective_permeability

  perm_dn_adj(:) = perm_dn

  Res = 0.d0

  v_darcy = 0.d0

  do iphase = 1,option%nphase - 1 ! No advection or diffusion in the salt phase

    ! Salt transport is either through aqueous or gas phases.
    if (iphase == PRECIPITATE_PHASE) cycle

    bc_type = ibndtype(iphase)
    ! Advection
    select case(bc_type)

    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC,DIRICHLET_SEEPAGE_BC)

      if (sco2_auxvar_up%mobility(iphase) + &
          sco2_auxvar_dn%mobility(iphase) > eps) then

      ! Harmonic mean on viscosity
      visc_mean = 2.d0 * (sco2_auxvar_up%visc(iphase) * &
                          sco2_auxvar_dn%visc(iphase)) / &
                         (sco2_auxvar_up%visc(iphase) + &
                          sco2_auxvar_dn%visc(iphase))

      ! STOMP takes harmonic mean on density when old velocity is 0

        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
        if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(SCO2_LIQUID_CONDUCTANCE_INDEX)
          case(GAS_PHASE)
            idof = auxvar_mapping(SCO2_GAS_CONDUCTANCE_INDEX)
          end select
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif

        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == SCO2_GAS_STATE) then
          boundary_pressure = sco2_auxvar_up%pres(option%gas_phase)
        else
          ! MAN: check if Pc is needed in the boundary auxvar
          boundary_pressure = sco2_auxvar_up%pres(iphase)

        endif

        density_kg_ave = SCO2AverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            sco2_auxvar_up%den_kg, &
                                            sco2_auxvar_dn%den_kg)

        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                         sco2_auxvar_dn%pres(iphase) + &
                         gravity_term

        if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
            ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              sco2_auxvar_up%pres(iphase) - &
              option%flow%reference_pressure < eps) then
            delta_pressure = 0.d0
          endif
        endif

        if (bc_type == DIRICHLET_SEEPAGE_BC) then
          if (delta_pressure < 0.d0) then
            delta_pressure = 0.d0
          endif
        endif

        dn_scale = 0.d0
        upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                                 PETSC_FALSE, &
                                 count_upwind_direction_flip_, &
                                 liq_upwind_flip_count_by_res, &
                                 liq_upwind_flip_count_by_jac)
        if (upwind) then
          mobility = sco2_auxvar_up%mobility(iphase)
          kr = sco2_auxvar_up%kr(iphase)
          xmass(:) = sco2_auxvar_up%xmass(:,iphase)
          uH = sco2_auxvar_up%H(iphase)
          density_kg_ave = sco2_auxvar_up%den_kg(iphase)
        else
          dn_scale = 1.d0
          mobility = sco2_auxvar_dn%mobility(iphase)
          kr = sco2_auxvar_dn%kr(iphase)
          xmass(:) = sco2_auxvar_dn%xmass(:,iphase)
          uH = sco2_auxvar_dn%H(iphase)
          density_kg_ave = sco2_auxvar_dn%den_kg(iphase)
        endif

        mobility = kr / visc_mean

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
      endif

    case(NEUMANN_BC)

      xmass(:) = 0.d0
      select case(iphase)
      case(LIQUID_PHASE)
        idof = auxvar_mapping(SCO2_LIQUID_FLUX_INDEX)
        if (ibndtype(SCO2_SALT_MASS_FRAC_DOF) == DIRICHLET_BC) then
          xmass(sid) = auxvars(SCO2_SALT_MASS_FRAC_DOF)
        else
          option%io_buffer = 'Salt concentration must be specified with a &
                              &DIRICHLET type BC.'
          call PrintErrMsg(option)
        endif
      case(GAS_PHASE)
        idof = auxvar_mapping(SCO2_GAS_FLUX_INDEX)
      end select

      ! MAN: might need to change how enthalpy is moved
      !      across the boundary, since this is a phase
      !      property and Neumann BC's are component-by-component
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then
          density_kg_ave = sco2_auxvar_up%den_kg(iphase)
          uH = sco2_auxvar_up%H(iphase)
        else
          dn_scale = 1.d0
          density_kg_ave = sco2_auxvar_dn%den_kg(iphase)
          uH = sco2_auxvar_dn%H(iphase)
        endif
      endif

    case default

      option%io_buffer = &
      'Boundary condition type not recognized in SCO2BCFlux phase loop.'
      call PrintErrMsg(option)

    end select

    if (dabs(v_darcy(iphase)) > 0.d0) then
      if (mobility > floweps ) then
        ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area
        ! mass_flux[kg phase/sec] = q[m^3 phase/sec] *
        !                             density_ave[kg phase/m^3 phase]
        tot_mass_flux = q * density_kg_ave
        ! comp_mass_flux[kg comp/sec] = tot_mass_flux[kg phase/sec] *
        !                                 xmass[kg comp/kg phase]
        do icomp = 1 , option%nflowspec - 1 ! Handle salt separately
          component_mass_flux = tot_mass_flux * xmass(icomp)
          Res(icomp) = Res(icomp) + component_mass_flux
        enddo
        if (.not. sco2_isothermal) then
          ! Energy flux
          Res(SCO2_ENERGY_EQUATION_INDEX) = Res(SCO2_ENERGY_EQUATION_INDEX) + &
                                          tot_mass_flux * uH
        endif
      endif
    endif


    ! Diffusion

    sat_dn = sco2_auxvar_dn%sat(iphase)
    if ((sat_dn > eps .and. ibndtype(iphase) /= NEUMANN_BC)) then
    ! Compute mole flux for CO2, mass flux for NaCl
    ! Harmonic diffusion coefficient
    den_up = sco2_auxvar_up%den(iphase)
    den_dn = sco2_auxvar_dn%den(iphase)
    den_kg_up = sco2_auxvar_up%den_kg(iphase)
    den_kg_dn = sco2_auxvar_dn%den_kg(iphase)

    dist_up = dist(0)
    dist_dn = dist(0)

    if (iphase == LIQUID_PHASE) then
      ! CO2 Mole Flux in the aqueous phase
      ! Include diffusion and longitudinal dispersion
      stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(co2_id,iphase) + &
                 sco2_auxvar_up%dispersivity(co2_id,iphase) * &
                 v_darcy(iphase))
      stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(co2_id,iphase) + &
                 sco2_auxvar_dn%dispersivity(co2_id,iphase) * &
                 v_darcy(iphase))

      ! Harmonic mean for diffusivity
      ! units = [kg/m^2/s bulk]
      stpd_ave_over_dist = stpd_up*stpd_dn / &
                           (5.d-1 * (stpd_up*dist_dn + stpd_dn*dist_up))

      ! units = kg/sec
      dtot_mole_flux_ddeltaX = stpd_ave_over_dist * area

      delta_xmol = sco2_auxvar_up%xmol(co2_id,iphase) * den_up - &
                   sco2_auxvar_dn%xmol(co2_id,iphase) * den_dn

      co2_mole_flux = dtot_mole_flux_ddeltaX * delta_xmol
    else
      ! Vapor Mole Flux in the gas phase
      ! Include diffusion and longitudinal dispersion
      stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(wid,iphase) + &
                 sco2_auxvar_up%dispersivity(wid,iphase) * &
                 v_darcy(iphase))
      stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(wid,iphase) + &
                 sco2_auxvar_dn%dispersivity(wid,iphase) * &
                 v_darcy(iphase))

      ! units = [kg/m^2/s bulk]
      stpd_ave_over_dist = stpd_up*stpd_dn / &
                           (5.d-1 * (stpd_up*dist_dn + stpd_dn*dist_up))

      ! units = kg/sec
      dtot_mole_flux_ddeltaX = stpd_ave_over_dist * area

      delta_xmol = sco2_auxvar_up%xmol(wid,iphase) * den_up - &
                   sco2_auxvar_dn%xmol(wid,iphase) * den_dn

      co2_mole_flux = -dtot_mole_flux_ddeltaX * delta_xmol
    endif

    ! Salt mass flux
    ! Patankar salt transport
    ! Include diffusion and longitudinal dispersion
    stpd_up = (sco2_auxvar_up%effective_diffusion_coeff(sid,iphase) + &
               sco2_auxvar_up%dispersivity(sid,iphase) * &
               v_darcy(iphase))
    stpd_dn = (sco2_auxvar_dn%effective_diffusion_coeff(sid,iphase) + &
               sco2_auxvar_dn%dispersivity(sid,iphase) * &
               v_darcy(iphase))

    tempreal = stpd_up*dist_up+stpd_dn*dist_dn
    if (tempreal > 0.d0) then
      stpd_ave_over_dist = stpd_up*stpd_dn / &
                           (5.d-1 * (stpd_up*dist_dn + stpd_dn*dist_up))
    else
      stpd_ave_over_dist = 0.d0
    endif

    al = max(-v_darcy(iphase),0.d0) + stpd_ave_over_dist * max((1.d0 - &
           (1.d-1 * dabs(v_darcy(iphase))/(stpd_ave_over_dist + epsilon))) ** 5, &
           0.d0)
    alp = max(v_darcy(iphase),0.d0) + stpd_ave_over_dist * max((1.d0 - &
            (1.d-1 * dabs(v_darcy(iphase))/(stpd_ave_over_dist + epsilon))) ** 5, &
            0.d0)

    salt_mass_flux = (alp * sco2_auxvar_up%xmass(sid,iphase) * &
                          sco2_auxvar_up%den_kg(iphase) - &
                          al * sco2_auxvar_dn%xmass(sid,iphase) * &
                          sco2_auxvar_dn%den_kg(iphase)) * area
    ! Diffusive component of salt flux
    ! units = kg/sec
    dsalt_mass_flux_ddeltaX = stpd_ave_over_dist * area

    delta_xmass =  sco2_auxvar_up%xmass(sid,iphase) * &
                   sco2_auxvar_up%den_kg(iphase) - &
                   sco2_auxvar_dn%xmass(sid,iphase) * &
                   sco2_auxvar_dn%den_kg(iphase)

    salt_diff_flux = dsalt_mass_flux_ddeltaX * delta_xmass

      if (iphase == LIQUID_PHASE) then
          water_mass_flux = -1.d0 * fmw_comp(1) * &
                   (co2_mole_flux + salt_diff_flux / fmw_comp(3))
      else
          water_mass_flux = -1.d0 * fmw_comp(1) * co2_mole_flux
      endif

      co2_mass_flux = fmw_comp(2) * co2_mole_flux

      Res(SCO2_WATER_EQUATION_INDEX) = Res(SCO2_WATER_EQUATION_INDEX) + &
                                       water_mass_flux
      Res(SCO2_CO2_EQUATION_INDEX) = Res(SCO2_CO2_EQUATION_INDEX) + &
                                     co2_mass_flux
      Res(SCO2_SALT_EQUATION_INDEX) = Res(SCO2_SALT_EQUATION_INDEX) + &
                                      salt_mass_flux
    endif
  enddo

  ! Conduction
  if (.not. sco2_isothermal) then
    heat_flux = 0.d0
    select case(ibndtype(SCO2_ENERGY_EQUATION_INDEX))
    case(DIRICHLET_BC)
      ! MAN: Need better thermal conductivity calculations, but right now
      !      taking roughly a weighted average of salt conductivity and
      !      pore/rock conductivity, assuming kgas ~ 0.
      sat_dn = sco2_auxvar_dn%sat(lid)

      ! derive wet and dry conductivities with anisotropy tensor and direction
      call thermal_cc_dn%thermal_conductivity_function% &
           TCondTensorToScalar(dist,option)
      call thermal_cc_dn%thermal_conductivity_function%CalculateTCond(sat_dn, &
             sco2_auxvar_dn%temp,sco2_auxvar_dn%effective_porosity, &
             k_eff_dn,dkeff_dn_dsatldn,dkeff_dn_dTdn,option)


      if (k_eff_dn > 0.d0) then
        k_eff_ave = k_eff_dn / dist(0)
      else
        k_eff_ave = 0.d0
      endif

      ! units:
      ! k_eff = W/K-m = J/s/K-m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = k_eff * delta_temp * area = J/s
      ! 1.0E-6 term accounts for change in units: J/s -> MJ/s

      delta_temp = sco2_auxvar_up%temp - sco2_auxvar_dn%temp
      dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
      heat_flux = dheat_flux_ddelta_temp * delta_temp
    case(NEUMANN_BC)
      ! Heat flux in MW/m^2
      heat_flux = auxvars(auxvar_mapping(SCO2_ENERGY_FLUX_INDEX)) * area
    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
          'SCO2BCFlux heat conduction loop.'
        call PrintErrMsg(option)
    end select

    ! ! MJ/s or MW
    Res(SCO2_ENERGY_EQUATION_INDEX) = Res(SCO2_ENERGY_EQUATION_INDEX) + heat_flux
  endif

end subroutine SCO2BCFlux

! ************************************************************************** !

subroutine SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                          sco2_auxvar_ss,sco2_auxvar,global_auxvar, &
                          global_auxvar_ss,material_auxvar, &
                          characteristic_curves, sco2_parameter, &
                          natural_id, scale,Res,aux_var_compute_only)
  !
  ! Computes source/sink terms for the residual
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use EOS_Water_module
  use EOS_Gas_module
  use Material_Aux_module
  use Characteristic_Curves_module

  implicit none

  type(option_type) :: option
  type(sco2_auxvar_type) :: sco2_auxvar,sco2_auxvar_ss
  type(global_auxvar_type) :: global_auxvar,global_auxvar_ss
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  type(sco2_parameter_type), pointer :: sco2_parameter
  PetscInt :: natural_id
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscBool :: aux_var_compute_only

  PetscReal :: qsrc(option%nflowdof)
  PetscInt :: flow_src_sink_type
  PetscReal :: mob_tot
  PetscReal :: xxss(option%nflowdof)
  PetscInt :: lid, gid, pid, wid, co2_id, co2_pressure_id, sid

  lid = option%liquid_phase
  gid = option%gas_phase
  pid = option%precipitate_phase
  co2_pressure_id = option%co2_pressure_id
  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  Res = 0.d0

  ! Extraction
  if (qsrc(wid)<0.d0 .or. qsrc(co2_id)<0.d0) then
    ! Use cell primary variables for state variable calculations
    xxss(SCO2_WATER_EQUATION_INDEX) = maxval(sco2_auxvar% &
                                    pres(option%liquid_phase:option%gas_phase))
    select case(global_auxvar%istate)
      case(SCO2_LIQUID_STATE)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar%xmass(co2_id,lid)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar%m_salt(1)
      case(SCO2_GAS_STATE)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar%pres(co2_pressure_id)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar%m_salt(2)
      case(SCO2_TRAPPED_GAS_STATE)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar%sat(gid)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar%m_salt(1)
      case(SCO2_LIQUID_GAS_STATE)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar%pres(gid)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar%m_salt(1)
    end select
    if (.not. sco2_isothermal) then
      xxss(SCO2_ENERGY_EQUATION_INDEX) = sco2_auxvar%temp
    endif
    global_auxvar_ss%istate = global_auxvar%istate
  else
    !Injection: use primary variables from user-supplied conditions
    select case(global_auxvar_ss%istate)
      case(SCO2_LIQUID_STATE)
        xxss(SCO2_WATER_EQUATION_INDEX) = sco2_auxvar_ss%pres(lid)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar_ss%xmass(co2_id,lid)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar_ss%m_salt(1)
      case(SCO2_GAS_STATE)
        xxss(SCO2_WATER_EQUATION_INDEX) = sco2_auxvar%pres(gid)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar_ss%pres(co2_pressure_id)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar_ss%m_salt(2)
      case(SCO2_TRAPPED_GAS_STATE)
        xxss(SCO2_WATER_EQUATION_INDEX) = sco2_auxvar%pres(lid)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar_ss%sat(gid)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar_ss%m_salt(1)
      case(SCO2_LIQUID_GAS_STATE)
        xxss(SCO2_WATER_EQUATION_INDEX) = sco2_auxvar%pres(lid)
        xxss(SCO2_CO2_EQUATION_INDEX) = sco2_auxvar_ss%pres(gid)
        xxss(SCO2_SALT_EQUATION_INDEX) = sco2_auxvar_ss%m_salt(1)
    end select
    if (.not. sco2_isothermal) then
      xxss(SCO2_ENERGY_EQUATION_INDEX) = sco2_auxvar_ss%temp
    endif
  endif

  call SCO2AuxVarCompute(xxss,sco2_auxvar_ss, global_auxvar_ss, &
                         material_auxvar, characteristic_curves, &
                         sco2_parameter, natural_id,option)

  if (aux_var_compute_only) return

  select case(flow_src_sink_type)
    case(TOTAL_MASS_RATE_SS)
      ! For extraction wells: apply a total mass sink and scale components
      ! by the mobility ratio. Stored in qsrc(1)
      mob_tot = sco2_auxvar%mobility(lid) + sco2_auxvar%mobility(gid)
      if (sco2_auxvar%sat(gid) <= 0.d0) then
        ! kg/sec total to kg/sec component
        Res(SCO2_WATER_EQUATION_INDEX) = qsrc(1) * sco2_auxvar%xmass(wid,lid)
        Res(SCO2_CO2_EQUATION_INDEX) = qsrc(1) * &
                                       sco2_auxvar%xmass(co2_id,lid)
        Res(SCO2_SALT_EQUATION_INDEX) = qsrc(1) * sco2_auxvar_ss%xmass(sid,lid)
        if (.not. sco2_isothermal) then
          Res(SCO2_ENERGY_EQUATION_INDEX) = qsrc(1) * sco2_auxvar_ss%H(lid)
        endif
      elseif (sco2_auxvar%sat(lid) <= 0.d0) then
        ! kg/sec total to kg/sec component
        Res(SCO2_WATER_EQUATION_INDEX) = qsrc(1) * sco2_auxvar%xmass(wid,gid)
        Res(SCO2_CO2_EQUATION_INDEX) = qsrc(1) * &
                                       sco2_auxvar%xmass(co2_id,gid)
        Res(SCO2_SALT_EQUATION_INDEX) = 0.d0
        if (.not. sco2_isothermal) then
          Res(SCO2_ENERGY_EQUATION_INDEX) = qsrc(1) * sco2_auxvar_ss%H(gid)
        endif
      else
        ! Water component
        Res(SCO2_WATER_EQUATION_INDEX) = qsrc(1) * &
                                         (sco2_auxvar%mobility(lid)/mob_tot * &
                                          sco2_auxvar%xmass(wid,lid) + &
                                          sco2_auxvar%mobility(gid)/mob_tot * &
                                          sco2_auxvar%xmass(wid,gid))
        ! CO2 component
        Res(SCO2_CO2_EQUATION_INDEX) = qsrc(1) * &
                                        (sco2_auxvar%mobility(lid)/mob_tot * &
                                         sco2_auxvar%xmass(co2_id,lid) + &
                                         sco2_auxvar%mobility(gid)/mob_tot * &
                                         sco2_auxvar%xmass(co2_id,gid))

        ! Salt component
        Res(SCO2_SALT_EQUATION_INDEX) = qsrc(1) * &
                                         (sco2_auxvar%mobility(lid)/mob_tot * &
                                          sco2_auxvar%xmass(sid,lid)) 
        if (.not. sco2_isothermal) then
          ! Energy                                  
          Res(SCO2_ENERGY_EQUATION_INDEX) = qsrc(1) * &
                                         (sco2_auxvar%mobility(lid)/mob_tot * &
                                          sco2_auxvar%H(lid) + &
                                          sco2_auxvar%mobility(gid)/mob_tot * &
                                          sco2_auxvar%H(gid))
        endif
        
      endif

    case(MASS_RATE_SS)
      Res(SCO2_WATER_EQUATION_INDEX) = qsrc(wid)
      Res(SCO2_CO2_EQUATION_INDEX) = qsrc(co2_id)
      Res(SCO2_SALT_EQUATION_INDEX) = qsrc(sid)
      if (.not. sco2_isothermal) then
        Res(SCO2_ENERGY_EQUATION_INDEX) = qsrc(wid) * &
                  sco2_auxvar_ss%H(lid) + qsrc(gid) * &
                  sco2_auxvar_ss%H(gid)
      endif
    case(SCALED_MASS_RATE_SS)
      Res(SCO2_WATER_EQUATION_INDEX) = qsrc(wid) * scale
      Res(SCO2_CO2_EQUATION_INDEX) = qsrc(co2_id) * scale
      Res(SCO2_SALT_EQUATION_INDEX) = qsrc(sid) * scale
      if (.not. sco2_isothermal) then
        Res(SCO2_ENERGY_EQUATION_INDEX) = scale * (qsrc(wid) * &
                  sco2_auxvar_ss%H(lid) + qsrc(gid) * &
                  sco2_auxvar_ss%H(gid))
      endif
    case(VOLUMETRIC_RATE_SS)
      ! This would have to be in m^3/sec phase
      Res(SCO2_WATER_EQUATION_INDEX) = qsrc(wid) * sco2_auxvar%den_kg(lid)
      Res(SCO2_CO2_EQUATION_INDEX) = qsrc(co2_id) * &
                                     sco2_auxvar%den_kg(gid)
      Res(SCO2_SALT_EQUATION_INDEX) = qsrc(sid) * SALT_DENSITY_KG
      if (.not. sco2_isothermal) then
        Res(SCO2_ENERGY_EQUATION_INDEX) = (qsrc(wid) * &
                   sco2_auxvar%den_kg(lid) * sco2_auxvar%H(lid) + qsrc(gid) * &
                   sco2_auxvar%den_kg(gid)) * sco2_auxvar%H(gid)
      endif
    case(SCALED_VOLUMETRIC_RATE_SS)
      ! This would have to be in m^3/sec phase
      Res(SCO2_WATER_EQUATION_INDEX) = qsrc(wid) * sco2_auxvar%den_kg(lid) * &
                                       scale
      Res(SCO2_CO2_EQUATION_INDEX) = qsrc(co2_id) * &
                                     sco2_auxvar%den_kg(gid) * scale
      
      Res(SCO2_SALT_EQUATION_INDEX) = qsrc(sid) * SALT_DENSITY_KG * &
                                      scale
      if (.not. sco2_isothermal) then
        Res(SCO2_ENERGY_EQUATION_INDEX) = (qsrc(wid) * &
                                        sco2_auxvar%H(lid) + qsrc(gid) * &
                                        sco2_auxvar%H(gid)) * scale
      endif
  end select

  ! If there's a heater
  if (.not. sco2_isothermal) Res(SCO2_ENERGY_EQUATION_INDEX) = &
                             Res(SCO2_ENERGY_EQUATION_INDEX) + qsrc(4)

end subroutine SCO2AuxVarComputeAndSrcSink

! ************************************************************************** !

subroutine SCO2AccumDerivative(sco2_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,J)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscInt :: idof, irow

  J = 0.d0
  res = 0.d0
  res_pert_plus = 0.d0
  res_pert_minus = 0.d0

  if (.not. sco2_central_diff_jacobian) then
    call SCO2Accumulation(sco2_auxvar(ZERO_INTEGER),global_auxvar, &
                           material_auxvar,soil_heat_capacity,option,res)
  endif


  if (sco2_central_diff_jacobian) then
    do idof = 1, option%nflowdof
      call SCO2Accumulation(sco2_auxvar(idof),global_auxvar, &
                            material_auxvar,soil_heat_capacity,option, &
                            res_pert_plus)

      call SCO2Accumulation(sco2_auxvar(idof+option%nflowdof), &
                            global_auxvar,material_auxvar,&
                            soil_heat_capacity,option, &
                            res_pert_minus)

       do irow = 1, option%nflowdof
         J(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0 * &
                         sco2_auxvar(idof)%pert)
       enddo !irow
    enddo ! idof
  else
    do idof = 1, option%nflowdof
      call SCO2Accumulation(sco2_auxvar(idof),global_auxvar, &
                            material_auxvar,soil_heat_capacity,option, &
                            res_pert_plus)

      do irow = 1, option%nflowdof
        J(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                        sco2_auxvar(idof)%pert
      enddo !irow
    enddo ! idof
  endif

end subroutine SCO2AccumDerivative

! ************************************************************************** !

subroutine SCO2FluxDerivative(sco2_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 thermal_cc_up, &
                                 sco2_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 thermal_cc_dn, &
                                 area, dist, upwind_direction_, &
                                 option,Jup,Jdn)
  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !
  use Option_module
  use Material_Aux_module
  use Upwind_Direction_module, only : count_upwind_direction_flip
  use Characteristic_Curves_Thermal_module

  implicit none

  type(sco2_auxvar_type) :: sco2_auxvar_up(0:), sco2_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  class(cc_thermal_type) :: thermal_cc_up, thermal_cc_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscInt :: idof, irow

  res = 0.d0
  res_pert_plus = 0.d0
  res_pert_minus = 0.d0

  Jup = 0.d0
  Jdn = 0.d0

  option%iflag = -2

  if (.not. sco2_central_diff_jacobian) then
    call SCO2Flux(sco2_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   thermal_cc_up, &
                   sco2_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   thermal_cc_dn, &
                   area,dist,upwind_direction_, &
                   option,v_darcy,res,&
                   PETSC_FALSE, & ! update the upwind direction
                   ! avoid double counting upwind direction flip
                   PETSC_FALSE) ! count upwind direction flip
  endif

  ! upgradient derivatives
  if(sco2_central_diff_jacobian) then
      do idof = 1, option%nflowdof
        call SCO2Flux(sco2_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up, &
                     thermal_cc_up, &
                     sco2_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_cc_dn, &
                     area,dist,upwind_direction_, &
                     option,v_darcy,res_pert_plus, &
                     PETSC_FALSE, & ! update the upwind direction
                     count_upwind_direction_flip)

        call SCO2Flux(sco2_auxvar_up(idof+option%nflowdof), &
                     global_auxvar_up,material_auxvar_up, &
                     thermal_cc_up, &
                     sco2_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_cc_dn, &
                     area,dist,upwind_direction_, &
                     option,v_darcy,res_pert_minus, &
                     PETSC_FALSE, & ! update the upwind direction
                     count_upwind_direction_flip)
        do irow = 1, option%nflowdof
          Jup(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/(2.d0 * &
                            sco2_auxvar_up(idof)%pert)
        enddo !irow
      enddo ! idof
  else
    do idof = 1, option%nflowdof
      call SCO2Flux(sco2_auxvar_up(idof),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_cc_up, &
                    sco2_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_cc_dn, &
                    area,dist,upwind_direction_, &
                    option,v_darcy,res_pert_plus, &
                    PETSC_FALSE, & ! update the upwind direction
                    count_upwind_direction_flip)

      do irow = 1, option%nflowdof
        Jup(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                         sco2_auxvar_up(idof)%pert
      enddo !irow
    enddo ! idof
  endif

  ! downgradient derivatives
  if (sco2_central_diff_jacobian) then
    do idof = 1, option%nflowdof
      call SCO2Flux(sco2_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_cc_up, &
                    sco2_auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_cc_dn, &
                    area,dist,upwind_direction_, &
                    option,v_darcy,res_pert_plus, &
                    PETSC_FALSE, & ! update the upwind direction
                    count_upwind_direction_flip)

      call SCO2Flux(sco2_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_cc_up, &
                    sco2_auxvar_dn(idof+option%nflowdof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_cc_dn, &
                    area,dist,upwind_direction_, &
                    option,v_darcy,res_pert_minus, &
                    PETSC_FALSE, & ! update the upwind direction
                    count_upwind_direction_flip)

      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0* &
                          sco2_auxvar_dn(idof)%pert)

      enddo !irow
    enddo ! idof
  else
    do idof = 1, option%nflowdof
      call SCO2Flux(sco2_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                    material_auxvar_up, &
                    thermal_cc_up, &
                    sco2_auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    thermal_cc_dn, &
                    area,dist,upwind_direction_, &
                    option,v_darcy,res_pert_plus, &
                    PETSC_FALSE, & ! update the upwind direction
                    count_upwind_direction_flip)

      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                          sco2_auxvar_dn(idof)%pert
      enddo !irow
    enddo ! idof
  endif
end subroutine SCO2FluxDerivative

! ************************************************************************** !

subroutine SCO2BCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                sco2_auxvar_up, &
                                global_auxvar_up, &
                                sco2_auxvar_dn,global_auxvar_dn, &
                                material_auxvar_dn, &
                                thermal_cc_dn, &
                                area,dist,upwind_direction_, &
                                option,Jdn)
  !
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use Material_Aux_module
  use Upwind_Direction_module, only : count_upwind_direction_flip
  use Characteristic_Curves_Thermal_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(SCO2_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(sco2_auxvar_type) :: sco2_auxvar_up, sco2_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  class(cc_thermal_type) :: thermal_cc_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscInt :: idof, irow

  res = 0.d0
  res_pert_plus = 0.d0
  res_pert_minus = 0.d0
  Jdn = 0.d0

  option%iflag = -2

  if (.not. sco2_central_diff_jacobian) then
    call SCO2BCFlux(ibndtype,auxvar_mapping,auxvars, &
                     sco2_auxvar_up,global_auxvar_up, &
                     sco2_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_cc_dn, &
                     area,dist,upwind_direction_, &
                     option,v_darcy,res, &
                     PETSC_FALSE, & ! update the upwind direction
                     ! avoid double counting upwind direction flip
                     PETSC_FALSE) ! count upwind direction flip
  endif


  ! downgradient derivatives
  if (sco2_central_diff_jacobian) then
    do idof = 1, option%nflowdof
      call SCO2BCFlux(ibndtype,auxvar_mapping,auxvars, &
                      sco2_auxvar_up,global_auxvar_up, &
                      sco2_auxvar_dn(idof),global_auxvar_dn, &
                      material_auxvar_dn, &
                      thermal_cc_dn, &
                      area,dist,upwind_direction_, &
                      option,v_darcy,res_pert_plus, &
                      PETSC_FALSE, & ! update the upwind direction
                      count_upwind_direction_flip)

      call SCO2BCFlux(ibndtype,auxvar_mapping,auxvars, &
                      sco2_auxvar_up,global_auxvar_up, &
                      sco2_auxvar_dn(idof+option%nflowdof),global_auxvar_dn, &
                      material_auxvar_dn, &
                      thermal_cc_dn, &
                      area,dist,upwind_direction_, &
                      option,v_darcy,res_pert_minus, &
                      PETSC_FALSE, & ! update the upwind direction
                      count_upwind_direction_flip)


      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0 * &
                          sco2_auxvar_dn(idof)%pert)
      enddo !irow
    enddo ! idof
  else
    do idof = 1, option%nflowdof
      call SCO2BCFlux(ibndtype,auxvar_mapping,auxvars, &
                      sco2_auxvar_up,global_auxvar_up, &
                      sco2_auxvar_dn(idof),global_auxvar_dn, &
                      material_auxvar_dn, &
                      thermal_cc_dn, &
                      area,dist,upwind_direction_, &
                      option,v_darcy,res_pert_plus, &
                      PETSC_FALSE, & ! update the upwind direction
                      count_upwind_direction_flip)

      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                          sco2_auxvar_dn(idof)%pert
      enddo !irow
    enddo ! idof
  endif

end subroutine SCO2BCFluxDerivative

! ************************************************************************** !

subroutine SCO2SrcSinkDerivative(option,source_sink,sco2_auxvar_ss, &
                                 sco2_auxvar,global_auxvar, &
                                 global_auxvar_ss,characteristic_curves, &
                                 sco2_parameter, natural_id,material_auxvar, &
                                 scale,Jac)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Option_module
  use Coupler_module
  use Characteristic_Curves_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  type(coupler_type), pointer :: source_sink
  type(sco2_auxvar_type) :: sco2_auxvar(0:), sco2_auxvar_ss(0:1)
  type(global_auxvar_type) :: global_auxvar, global_auxvar_ss
  class(characteristic_curves_type) :: characteristic_curves
  type(sco2_parameter_type), pointer :: sco2_parameter
  PetscInt :: natural_id
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)

  PetscReal :: qsrc(option%nflowdof)
  PetscInt :: flow_src_sink_type
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscInt :: idof, irow

  res = 0.d0
  res_pert_plus = 0.d0
  res_pert_minus = 0.d0
  Jac = 0.d0

  qsrc = source_sink%flow_condition%sco2%rate%dataset%rarray(:)
  flow_src_sink_type = source_sink%flow_condition%sco2%rate%itype

  option%iflag = -3

  if (.not. sco2_central_diff_jacobian) then
    ! Index 0 contains user-specified conditions
    ! Index 1 contains auxvars to be used in src/sink calculations
    call SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                    sco2_auxvar_ss(ZERO_INTEGER), sco2_auxvar(ZERO_INTEGER), &
                    global_auxvar, global_auxvar_ss, material_auxvar, &
                    characteristic_curves, sco2_parameter, natural_id, &
                    scale,res,PETSC_FALSE)
  endif


  ! downgradient derivatives
  if (sco2_central_diff_jacobian) then
    do idof = 1, option%nflowdof
      call SCO2AuxVarCopy(sco2_auxvar_ss(ZERO_INTEGER), &
                             sco2_auxvar_ss(ONE_INTEGER), option)
      call SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                        sco2_auxvar_ss(ONE_INTEGER), &
                        sco2_auxvar(idof), global_auxvar, global_auxvar_ss, &
                        material_auxvar,characteristic_curves, &
                        sco2_parameter, natural_id, scale, res_pert_plus, &
                        PETSC_FALSE)

      call SCO2AuxVarCopy(sco2_auxvar_ss(ZERO_INTEGER), &
                             sco2_auxvar_ss(ONE_INTEGER), option)
      call SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                        sco2_auxvar_ss(ONE_INTEGER), &
                        sco2_auxvar(idof+option%nflowdof),global_auxvar,&
                        global_auxvar_ss, &
                        material_auxvar,characteristic_curves, &
                        sco2_parameter, natural_id, scale, res_pert_minus, &
                        PETSC_FALSE)

      do irow = 1, option%nflowdof
        Jac(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0 * &
                          sco2_auxvar(idof)%pert)
      enddo !irow
    enddo ! idof
  else
    do idof = 1, option%nflowdof
      call SCO2AuxVarCopy(sco2_auxvar_ss(ZERO_INTEGER), &
                             sco2_auxvar_ss(ONE_INTEGER), option)
      call SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                        sco2_auxvar_ss(ONE_INTEGER), &
                        sco2_auxvar(idof),global_auxvar, global_auxvar_ss, &
                        material_auxvar,characteristic_curves, &
                        sco2_parameter, natural_id, scale, res_pert_plus, &
                        PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jac(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                          sco2_auxvar(idof)%pert
      enddo !irow
    enddo ! idof
  endif

end subroutine SCO2SrcSinkDerivative

! ************************************************************************** !

function SCO2AverageDensity(iphase,istate_up,istate_dn,density_up,density_dn)
  !
  ! Averages density, using opposite cell density if phase non-existent
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: SCO2AverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0

  if (iphase == LIQUID_PHASE) then
    if (istate_up == SCO2_GAS_STATE) then
      SCO2AverageDensity = density_dn(iphase)
    else if (istate_dn == SCO2_GAS_STATE) then
      SCO2AverageDensity = density_up(iphase)
    else
      SCO2AverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == SCO2_LIQUID_STATE) then
      SCO2AverageDensity = density_dn(iphase)
    else if (istate_dn == SCO2_LIQUID_STATE) then
      SCO2AverageDensity = density_up(iphase)
    else
      SCO2AverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
    endif
  endif

end function SCO2AverageDensity

! ************************************************************************** !

end module SCO2_Common_module
