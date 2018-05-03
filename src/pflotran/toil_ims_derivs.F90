module TOilIms_derivs_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use AuxVars_TOilIms_module 

  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use PM_TOilIms_Aux_module

  implicit none
  
  private 
#define TOIL_CONVECTION
#define TOIL_CONDUCTION

! Cutoff parameters - no public
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  public :: toil_accum_derivs_alyt, &
            TOilImsFluxPFL_derivs, &
            TOilImsAverageDensity_derivs, &
            MoleFluxDerivs, &
            EnergyDrivenFluxDerivs, &
            DeltaPressureDerivs_up_and_down, &
            V_Darcy_Derivs, &
            InjectionEnergyPartDerivs, &
            Qsrc_mol_derivs

contains


! ************************************************************************** !

!subroutine InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, &
                               !r, enthalpy, &
                               !hw_dp, hw_dT)

subroutine InjectionEnergyPartDerivs(d_inj_en_part, denth_bool, r, &
                                      enthalpy, hw_dp, hw_dT)

  implicit none
  PetscReal, dimension(1:3) :: d_inj_en_part
  PetscReal :: denth_bool, r, enthalpy, hw_dp, hw_dT
     
  !r *  enthalpy

  d_inj_en_part = 0.d0

  !! w.r.t. oil pressure
  d_inj_en_part(1) = denth_bool * r * hw_dp


  !! w.r.t. temperature
  d_inj_en_part(3) = denth_bool * r * hw_dT

 

end subroutine InjectionEnergyPartDerivs

! ************************************************************************** !

subroutine Qsrc_mol_derivs(d_qsrc_mol,dden_bool, qsrc, dden_dp, dden_dt, sc)
  
  implicit none

  PetscReal, dimension(1:3) :: d_qsrc_mol
  PetscReal :: dden_bool
  PetscReal :: qsrc, dden_dp, dden_dt, sc

  !! qsrc_mol = sc*qsrc*den

  d_qsrc_mol = 0.d0

  !! w.r.t. pressure:
  d_qsrc_mol(1) = dden_bool * qsrc * dden_dp

  !! w.r.t. saturation
  !!  (is 0)

  !! w.r.t. temperature:
  d_qsrc_mol(3) = dden_bool * qsrc * dden_dT

  d_qsrc_mol = d_qsrc_mol * sc


end subroutine Qsrc_mol_derivs

! ************************************************************************** !

subroutine TOilImsFluxPFL_derivs(toil_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       sir_up, &
                       thermal_conductivity_up, &
                       toil_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area, dist, &
                       option,v_darcy,Res, &
                       jup, jdn)

  use Option_module
  use Material_Aux_class
  use Connection_module
 
  ! no fractures considered for now
  ! use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn
  PetscReal, dimension(1:3,1:3) :: jup, jdn
  class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  !type(toil_ims_parameter_type) :: parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  !PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscInt :: energy_id
  PetscInt :: iphase
 
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn           ! no mole fractions
  PetscReal :: delta_pressure, delta_temp !, delta_xmol,

  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux

  ! no diff fluxes - arrays used for debugging only
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  
  PetscReal :: dummy_perm_up, dummy_perm_dn

  PetscReal :: up_scale, dn_scale

  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: d_delta_pres_dp_up, d_delta_pres_dp_dn
  PetscReal :: d_delta_pres_dT_up, d_delta_pres_dT_dn

  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: d_delta_temp_dt_up, d_delta_temp_dt_dn, dheat_flux_ddelta_temp
  PetscReal :: d_delta_pres_ds_up, d_delta_pres_ds_dn   

  PetscReal, dimension(1:3) :: d_v_darcy_up, d_v_darcy_dn
  PetscReal, dimension(1:3) :: d_q_up, d_q_dn
  PetscReal, dimension(1:3) :: d_mole_flux_up, d_mole_flux_dn
  PetscReal, dimension(1:3) :: d_energy_flux_up, d_energy_flux_dn



  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
      
  Res = 0.d0
  
  v_darcy = 0.d0

  jup = 0.d0
  jdn = 0.d0


#ifdef TOIL_CONVECTION
  do iphase = 1, option%nphase
    !print *, "iphase: ", iphase, ", phases: ", option%nphase
    d_v_darcy_up = 0.d0
    d_v_darcy_dn = 0.d0
    d_q_up = 0.d0
    d_q_dn = 0.d0
    d_mole_flux_up = 0.d0
    d_mole_flux_dn = 0.d0
    d_energy_flux_up = 0.d0
    d_energy_flux_dn = 0.d0
 
    if (toil_auxvar_up%mobility(iphase) + &
        toil_auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    ! an alternative could be to avergae using oil_sat
    !density_kg_ave = 0.5d0* ( toil_auxvar_up%den_kg(iphase) + &
    !                          toil_auxvar_dn%den_kg(iphase) )
    density_kg_ave = TOilImsAverageDensity_derivs(toil_auxvar_up%sat(iphase), &
                     toil_auxvar_dn%sat(iphase), &
                     toil_auxvar_up%den_kg(iphase), &
                     toil_auxvar_dn%den_kg(iphase), &
                     ddensity_kg_ave_dden_kg_up, &
                     ddensity_kg_ave_dden_kg_dn)


    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = toil_auxvar_up%pres(iphase) - &
                     toil_auxvar_dn%pres(iphase) + &
                     gravity_term



!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_dphi(iphase) = delta_pressure
!#endif

    ! upwinding the mobilities and enthalpies
    up_scale = 0.d0 !! ADDED
    dn_scale = 0.d0 !! ADDED
    if (delta_pressure >= 0.D0) then
      up_scale = 1.d0
      mobility = toil_auxvar_up%mobility(iphase)
      H_ave = toil_auxvar_up%H(iphase)
      uH = H_ave
#ifdef TOIL_DEN_UPWIND
      density_ave = toil_auxvar_up%den(iphase)
      ddensity_ave_dden_up = 1.d0
      ddensity_ave_dden_dn = 0.d0
#endif
    else
      dn_scale = 1.d0
      mobility = toil_auxvar_dn%mobility(iphase)
      H_ave = toil_auxvar_dn%H(iphase)
      uH = H_ave
#ifdef TOIL_DEN_UPWIND
      density_ave = toil_auxvar_dn%den(iphase)
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 1.d0
#endif
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure



      ! if comments below, use upwinding value
      !density_ave = 0.5d0*( toil_auxvar_up%den(iphase) + &
      !                      toil_auxvar_dn%den(iphase))
#ifndef TOIL_DEN_UPWIND
      density_ave = TOilImsAverageDensity_derivs(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den(iphase), &
                           toil_auxvar_dn%den(iphase), &
                           ddensity_ave_dden_up, &
                           ddensity_ave_dden_dn)
#endif 

      !! defer delta pressure derivtives to here because we know density average 
      !! derivatives will have been calculated by this point
      call  DeltaPressureDerivs_up_and_down(d_delta_pres_dp_up, d_delta_pres_dp_dn,     &
                                            d_delta_pres_dT_up, d_delta_pres_dT_dn,     &
                                            dist_gravity,                               &
                                            ddensity_ave_dden_up, ddensity_ave_dden_dn, &
                                            toil_auxvar_up%d%dden_dp(iphase,1),                &
                                            toil_auxvar_dn%d%dden_dp(iphase,1),                &
                                            toil_auxvar_up%d%dden_dt(iphase),                 &
                                            toil_auxvar_dn%d%dden_dt(iphase),                 &
                                            toil_auxvar_up%d%dp_dsat(iphase),                 &
                                            toil_auxvar_dn%d%dp_dsat(iphase),                 &
                                            d_delta_pres_ds_up, d_delta_pres_ds_dn,&
                                            toil_ims_fmw_comp(iphase))

                                           !dp_ds_up, dp_ds_dn,                         &
                                           !ddelta_pressure_dsatup, ddelta_pressure_dsatdn )


      call v_darcy_derivs(d_v_darcy_up, toil_auxvar_up%d%dmobility(iphase, 1:3),     &
                          up_scale, delta_pressure, mobility, perm_ave_over_dist(iphase),  &
                          d_delta_pres_dp_up, d_delta_pres_dt_up, d_delta_pres_ds_up)
      call v_darcy_derivs(d_v_darcy_dn, toil_auxvar_dn%d%dmobility(iphase, 1:3),     &
                          dn_scale, delta_pressure, mobility, perm_ave_over_dist(iphase),  &
                          d_delta_pres_dp_dn, d_delta_pres_dt_dn, d_delta_pres_ds_dn)

      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  

      d_q_up = d_v_darcy_up * area
      d_q_dn = d_v_darcy_dn * area

      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave

      call MoleFluxDerivs(d_mole_flux_up, d_q_up, q, density_ave, ddensity_ave_dden_up, &
                          toil_auxvar_up%d%dden_dp(iphase,1), toil_auxvar_up%d%dden_dt(iphase))
      call MoleFluxDerivs(d_mole_flux_dn, d_q_dn, q, density_ave, ddensity_ave_dden_dn, &
                          toil_auxvar_dn%d%dden_dp(iphase,1), toil_auxvar_dn%d%dden_dt(iphase))

      ! Res[kmol total/sec]

      ! Res[kmol phase/sec] = mole_flux[kmol phase/sec]  
      Res(iphase) = Res(iphase) + mole_flux 
      jup(iphase,1:3) = jup(iphase, 1:3) + d_mole_flux_up
      jdn(iphase, 1:3) = jdn(iphase, 1:3) + d_mole_flux_dn

      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/kmol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo

!#ifdef DEBUG_FLUXES  
!      do icomp = 1, option%nflowspec
!        adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
!      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      do icomp = 1, option%nflowspec
!        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
!#endif

      Res(energy_id) = Res(energy_id) + mole_flux * uH


      !! by `energy driven flux' mean the term moleFlux * uH
      call EnergyDrivenFluxDerivs(d_energy_flux_up, d_mole_flux_up, uH, up_scale, mole_flux, &
                            toil_auxvar_up%d%dH_dp(iphase), toil_auxvar_up%d%dH_dt(iphase))
      call EnergyDrivenFluxDerivs(d_energy_flux_dn, d_mole_flux_dn, uH, dn_scale, mole_flux, &
                            toil_auxvar_dn%d%dH_dp(iphase), toil_auxvar_dn%d%dH_dt(iphase))
      jup(energy_id, 1:3) = jup(energy_id, 1:3) + d_energy_flux_up
      jdn(energy_id, 1:3) = jdn(energy_id, 1:3) + d_energy_flux_dn


!#ifdef DEBUG_FLUXES  
!      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_dphi(iphase) = delta_pressure
!      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
!#endif

    endif  ! if mobility larger than given tolerance                 

  enddo
#endif 
! TOIL_CONVECTION

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then  
!    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
!    write(debug_unit,'(a,7es24.15)') 'adv flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'adv flux (gas):', debug_flux(:,2)
!  endif
!  debug_flux = 0.d0
!#endif                    

#ifdef TOIL_CONDUCTION
  ! model for liquid + gas
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
  !k_eff_up = thermal_conductivity_up(1) + &
  !           sqrt(gen_auxvar_up%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  !k_eff_dn = thermal_conductivity_dn(1) + &
  !           sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  !if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
  !  k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  !else
  !  k_eff_ave = 0.d0
  !endif
  ! considered the formation fully saturated in water for heat conduction 
  k_eff_up = thermal_conductivity_up(1)
  k_eff_dn = thermal_conductivity_dn(1)
  if (k_eff_up > 0.d0 .or. k_eff_dn > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif

  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = toil_auxvar_up%temp - toil_auxvar_dn%temp

  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s

  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux

  !!  analytical derivatives:
  d_delta_temp_dt_up = 1.d0
  d_delta_temp_dt_dn = - 1.d0

  dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s

  jup(energy_id, 3) = jup(energy_id, 3) + d_delta_temp_dt_up*dheat_flux_ddelta_temp
  jdn(energy_id, 3) = jdn(energy_id, 3) + d_delta_temp_dt_dn*dheat_flux_ddelta_temp


! CONDUCTION
#endif

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
!  if (debug_flag > 0) then  
!    write(debug_unit,'(a,7es24.15)') 'dif flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'dif flux (gas):', debug_flux(:,2)
!  endif
!#endif

end subroutine TOilImsFluxPFL_derivs

! ************************************************************************** !

function TOilImsAverageDensity_derivs(sat_up,sat_dn,density_up,density_dn, dden_up, dden_dn)
  ! 
  ! Modification of TOilImsAverageDensity which computes derivatives
  ! Daniel Stone, March-May 2018
  ! 

  implicit none

  PetscReal :: sat_up, sat_dn
  PetscReal :: density_up, density_dn, dden_up, dden_dn

  PetscReal :: TOilImsAverageDensity_derivs

  dden_up = 0.d0
  dden_dn = 0.d0

  if (sat_up < eps ) then
    TOilImsAverageDensity_derivs = density_dn
    dden_dn = 1.d0
  else if (sat_dn < eps ) then 
    TOilImsAverageDensity_derivs = density_up
    dden_up = 1.d0
  else ! in here we could use an armonic average, 
       ! other idea sat weighted average but it needs truncation
    TOilImsAverageDensity_derivs = 0.5d0*(density_up+density_dn)
    dden_up = 0.5d0
    dden_dn = 0.5d0
  end if

end function TOilImsAverageDensity_derivs

! ************************************************************************** !

subroutine DeltaPressureDerivs_up_and_down(ddelta_pressure_dpup, ddelta_pressure_dpdn, &
                                           ddelta_pressure_dTup, ddelta_pressure_dTdn, &
                                           dist_gravity,                               &
                                           ddensity_ave_dden_up, ddensity_ave_dden_dn, &
                                           dden_dp_up, dden_dp_dn,                     &
                                           dden_dt_up, dden_dt_dn,                     &
                                           dp_ds_up, dp_ds_dn,                         &
                                           ddelta_pressure_dsatup, ddelta_pressure_dsatdn, &
                                           fmw)


  implicit none

  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  PetscReal :: ddelta_pressure_dsatup, ddelta_pressure_dsatdn
  PetscReal :: dist_gravity
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dden_dp_up, dden_dp_dn
  PetscReal :: dden_dt_up, dden_dt_dn
  PetscReal :: dp_ds_up, dp_ds_dn
  PetscReal :: fmw
  PetscReal :: fmw_use


  !fmw_use = 0.d0
  fmw_use = fmw


  !gravity_term = density_kg_ave * dist_gravity
  !delta_pressure = toil_auxvar_up%pres(iphase) - &
                   !toil_auxvar_dn%pres(iphase) + &
                   !gravity_term

  ! w.r.t. pressure:
  ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                         ddensity_ave_dden_up * &
                         dden_dp_up * fmw_use
  ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                         ddensity_ave_dden_dn * &
                         dden_dp_up * fmw_use

  ! w.r.t. saturation:
  ddelta_pressure_dsatup = dp_ds_up

  ddelta_pressure_dsatdn = -1.0*dp_ds_dn

  ! w.r.t. temperature:
  ddelta_pressure_dTup = dist_gravity * &
                         ddensity_ave_dden_up * &
                         dden_dt_up * fmw_use
  ddelta_pressure_dTdn = dist_gravity * &
                         ddensity_ave_dden_dn * &
                         dden_dt_dn * fmw_use


end subroutine DeltaPressureDerivs_up_and_down

! ************************************************************************** !

subroutine EnergyDrivenFluxDerivs(d_energy_flux, d_mole_flux, uH, updn_scale, mole_flux, &
                            dH_dp, dH_dt)

  implicit none

  PetscReal, dimension(1:3) :: d_energy_flux
  PetscReal, dimension(1:3) :: d_mole_flux
  PetscReal :: uH, updn_scale, mole_flux
  PetscReal :: dH_dp, dH_dt

  !Res(energy_id) = Res(energy_id) + mole_flux * uH

  ! w.r.t. oil pressure
  d_energy_flux(1) = d_mole_flux(1)*uH + updn_scale*mole_flux*dH_dp

  ! w.r.t. oil saturation
  d_energy_flux(2) = d_mole_flux(2)*uH

  ! w.r.t. temperature
  d_energy_flux(3) = d_mole_flux(3)*uH + updn_scale*mole_flux*dH_dt

end subroutine EnergyDrivenFluxDerivs

! ************************************************************************** !

subroutine v_darcy_derivs(d_v_darcy, dmobility, updn_scale, delta_pressure, mobility, &
                          perm_ave_over_dist, d_delta_pres_dp, d_delta_pres_dt, &
                          d_delta_pres_ds)

  implicit none

  PetscReal, dimension(1:3) ::  d_v_darcy
  PetscReal, dimension(1:3) ::  dmobility
  PetscReal :: updn_scale, delta_pressure, mobility, perm_ave_over_dist
  PetscReal :: d_delta_pres_dp, d_delta_pres_dt
  PetscReal :: d_delta_pres_ds


  !v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure

   d_v_darcy = 0.d0

  ! w.r.t. oil pressure
  d_v_darcy(1) = updn_scale*dmobility(1)*delta_pressure + mobility*d_delta_pres_dp

  ! w.r.t. oil saturation
  d_v_darcy(2) = updn_scale*dmobility(2)*delta_pressure + mobility*d_delta_pres_ds

  ! w.r.t. temperature 
  d_v_darcy(3) = updn_scale*dmobility(3)*delta_pressure + mobility*d_delta_pres_dt


  !! scale by perm
  d_v_darcy = perm_ave_over_dist*d_v_darcy

end subroutine V_Darcy_Derivs

! ************************************************************************** !

subroutine MoleFluxDerivs(d_mole_flux, d_q, q, density_ave, ddensity_ave_dden, &
                          dden_dp, dden_dt)

  implicit none
  
  PetscReal, dimension(1:3) :: d_mole_flux
  PetscReal, dimension(1:3) :: d_q
  PetscReal :: q, density_ave
  PetscReal :: ddensity_ave_dden, dden_dp, dden_dt


  !mole_flux = q*density_ave
  ! d_mole_flux(i) = deriv of mole flux w.r.t. variable i
  ! 
  ! variables are:
  !
  ! pres
  ! sat 
  ! temp

  ! assume access to 
  !
  ! q, density ave
  ! dq: derivative array
  ! 
  ! ddensity_ave_dden
  ! so can chain rule with other density derivatives

  d_mole_flux = 0.d0

  ! w.r.t. oil pressure
  d_mole_flux(1) =  d_q(1)*density_ave + q*ddensity_ave_dden*dden_dp
  !                                       (  d (ave den) / d (p)    )

  ! w.r.t. oil sat
  d_mole_flux(2) = d_q(2)*density_ave

  ! w.r.t. temperature 
  d_mole_flux(3) =  d_q(3)*density_ave + q*ddensity_ave_dden*dden_dt
  !                                       (  d (ave den) / d (t)    )

end subroutine MoleFluxDerivs

! ************************************************************************** !

subroutine toil_accum_derivs_alyt(toil_auxvar,material_auxvar, option, j, soil_heat_capacity)

  use Material_Aux_class
  use Option_module

  implicit none

  !! Inputs:
  type(auxvar_toil_ims_type) :: toil_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  !! Outputs:
  PetscReal, dimension(1:3,1:3) :: j !! entry j is partial accum(j) / partial v_j
  !! workers:
  PetscReal :: porosity, volume
  PetscInt :: oid, lid, energy_id, iphase

  oid = option%oil_phase
  lid = option%liquid_phase
  energy_id = option%energy_id
  volume = material_auxvar%volume
  porosity = toil_auxvar%effective_porosity


  j = 0.d0


  !!    OIL EQUATION:

  !! w.r.t. pressure
  !print *, "ddendp: ", toil_auxvar%d%dden_dp(oid,1)
  j(oid, 1) = porosity*toil_auxvar%sat(oid)*toil_auxvar%d%dden_dp(oid,1) + &
              toil_auxvar%d%dpor_dp*(toil_auxvar%sat(oid)*toil_auxvar%den(oid))

  !! w.r.t. sat:
  j(oid, 2)  = porosity*toil_auxvar%den(oid)

  !! w.r.t. temp:
  j(oid,3) = porosity*toil_auxvar%sat(oid)*toil_auxvar%d%dden_dT(oid)


  !! END OIL EQUATION


  !!     LIQUID EQUATION

  !! w.r.t. pressure:
  !print *, "ddendp: ", toil_auxvar%d%dden_dp(lid,1)
  j(lid,1)  = porosity*toil_auxvar%sat(lid)*toil_auxvar%d%dden_dp(lid,1) + &
              toil_auxvar%d%dpor_dp*(toil_auxvar%sat(lid)*toil_auxvar%den(lid))

  !! w.r.t. sat:
  j(lid,2)  = -1.d0*toil_auxvar%den(lid)*porosity

  !! w.r.t. temp
  j(lid,3) =  porosity*toil_auxvar%sat(lid)*toil_auxvar%d%dden_dT(lid)

  !! END LIQUID EQUATION

  !! first a sum over the two phases, with the term being
  !! 
  !!    (poro) (sat) (den) (U)

  !!  liquid phase
  !! 
  !! w.r.t pressure
  j(energy_id, 1) = j(energy_id, 1) + & 
                    porosity * &
                    toil_auxvar%sat(lid)* ( &
                    toil_auxvar%d%dden_dp(lid,1)*toil_auxvar%U(lid) + &
                    toil_auxvar%den(lid)*toil_auxvar%d%dU_dp(lid) )
  !! and density w.r.t. pressure:
  j(energy_id, 1) = j(energy_id, 1) + & 
                    toil_auxvar%d%dpor_dp*toil_auxvar%sat(lid)*toil_auxvar%den(lid)*toil_auxvar%u(lid)

  !! w.r.t oil sat:
  j(energy_id, 2) = j(energy_id, 2) - & !! note negative, next term is scaled by dsl/dso 
                    porosity*toil_auxvar%den(lid)*toil_auxvar%U(lid)

  !! w.r.t. temp
  j(energy_id,3) = j(energy_id,3) + &
                   porosity * &
                   toil_auxvar%sat(lid)* ( &
                   toil_auxvar%d%dden_dt(lid)*toil_auxvar%U(lid) + &
                   toil_auxvar%den(lid)*toil_auxvar%d%dU_dT(lid)  )

  !!  oil phase
  !!
  !! w.r.t pressure
  j(energy_id, 1) = j(energy_id, 1) + & 
                    porosity * &
                    toil_auxvar%sat(oid)* ( &
                    toil_auxvar%d%dden_dp(oid,1)*toil_auxvar%U(oid) + &
                    toil_auxvar%den(oid)*toil_auxvar%d%dU_dp(oid) )
  !! and density w.r.t. pressure:
  j(energy_id, 1) = j(energy_id, 1) + & 
                    toil_auxvar%d%dpor_dp*toil_auxvar%sat(oid)*toil_auxvar%den(oid)*toil_auxvar%u(oid)

  !! w.r.t oil sat:
  j(energy_id, 2) = j(energy_id, 2) + & 
                    porosity*toil_auxvar%den(oid)*toil_auxvar%U(oid)

  !! w.r.t. temp
  j(energy_id,3) = j(energy_id,3) + &
                   porosity * &
                   toil_auxvar%sat(oid)* ( &
                   toil_auxvar%d%dden_dt(oid)*toil_auxvar%U(oid) + &
                   toil_auxvar%den(oid)*toil_auxvar%d%dU_dT(oid)  )

  !! also the (1-por) ... term
  j(energy_id,1) = j(energy_id,1) - toil_auxvar%d%dpor_dp*material_auxvar%soil_particle_density*soil_heat_capacity*toil_auxvar%temp
  j(energy_id,3) = j(energy_id,3) + (1.d0 - porosity)*material_auxvar%soil_particle_density * &
                                                             soil_heat_capacity

  !! END ENERGY EQUATION


j = j*volume

end subroutine toil_accum_derivs_alyt

! ************************************************************************** !
end module TOilIms_derivs_module
