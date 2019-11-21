module Hydrate_Common_module

  use Hydrate_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use petscsys

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

#define CONVECTION
#define LIQUID_DARCY_FLUX
#define GAS_DARCY_FLUX
#define DIFFUSION
#define LIQUID_DIFFUSION
#define GAS_DIFFUSION
#define CONDUCTION

#define WATER_SRCSINK
#define AIR_SRCSINK
#define ENERGY_SRCSINK
  

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24


  public :: HydrateAccumulation, &
            HydrateFlux, &
            HydrateBCFlux, &
            HydrateSrcSink, &
            HydrateAccumDerivative, &
            HydrateFluxDerivative, &
            HydrateBCFluxDerivative, &
            HydrateSrcSinkDerivative
            
  public :: HydrateDiffJacobian, &
            HydrateAuxVarDiff

contains

! ************************************************************************** !

subroutine HydrateAccumulation(hyd_auxvar,global_auxvar,material_auxvar, &
                               z,offset,methanogenesis,soil_heat_capacity, &
                               option,Res,Jac,analytical_derivatives,debug_cell)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual, for the hydrate sub-pm
  !
  ! Author: Michael Nole
  ! Date: 03/01/19
  !

  use Option_module
  use Material_Aux_class

  implicit none

  type(hydrate_auxvar_type) :: hyd_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: z, offset
  type(methanogenesis_type), pointer :: methanogenesis
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: debug_cell

  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase

  PetscReal :: porosity, volume
  PetscReal :: volume_over_dt
  PetscReal :: q_meth

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  volume_over_dt = material_auxvar%volume / option%flow_dt
  porosity = hyd_auxvar%effective_porosity
  volume = material_auxvar%volume
  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] *
    !                           den[kmol phase/m^3 phase] *
    !                           xmol[kmol comp/kmol phase]
    do icomp = 1, option%nflowspec
      Res(icomp) = Res(icomp) + hyd_auxvar%sat(iphase) * &
                                hyd_auxvar%den(iphase) * &
                                hyd_auxvar%xmol(icomp,iphase)
    enddo
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] *
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
    Res(energy_id) = Res(energy_id) + hyd_auxvar%sat(iphase) * &
                                      hyd_auxvar%den(iphase) * &
                                      hyd_auxvar%U(iphase)
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] +
  !                (1-por)[m^3 rock/m^3 bulk] *
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * hyd_auxvar%temp) * volume_over_dt

  q_meth = 0.d0
  if (associated(methanogenesis) .and. offset > 0.d0) then
    if (material_auxvar%id /= 1000) then
      call HydrateMethanogenesis(z, offset, methanogenesis, q_meth)
      !kmol/m^3/s to kmol/s
      q_meth = q_meth*(1.d0 - porosity)*volume
      Res(air_comp_id) = Res(air_comp_id) + q_meth
    endif
  endif
end subroutine HydrateAccumulation


! ************************************************************************** !

subroutine HydrateFlux(hyd_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area, dist, upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,Res,Jup,Jdn, &
                       analytical_derivatives, &
                       update_upwind_direction_, &
                       count_upwind_direction_flip_, &
                       debug_connection)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Fracture_module
  use Klinkenberg_module
  use Upwind_Direction_module
  
  implicit none
  
  type(hydrate_auxvar_type) :: hyd_auxvar_up, hyd_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(methanogenesis_type) :: methanogenesis
  type(hydrate_parameter_type) :: hydrate_parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: iphase
  
  PetscReal :: xmol(option%nflowspec)
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass
  PetscReal :: delta_X_whatever
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, q
  PetscReal :: tot_mole_flux, wat_mole_flux, air_mole_flux
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  
  PetscReal :: dummy_dperm_up, dummy_dperm_dn
  PetscReal :: temp_perm_up, temp_perm_dn

  ! Darcy flux
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dpaup, ddelta_pressure_dpadn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  
  PetscReal :: up_scale, dn_scale
  PetscBool :: upwind
  PetscReal :: tot_mole_flux_ddel_pressure
  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dtot_mole_flux_dp, dtot_mole_flux_dT, dtot_mole_flux_dsatg
  PetscReal :: dpl_dsatg
  PetscReal :: ddelta_pressure_pl
  PetscReal :: perm3(3)

  ! Diffusion
  PetscReal :: dstpd_up_dporup, dstpd_dn_dpordn
  PetscReal :: dstpd_up_dsatup, dstpd_dn_dsatdn
  PetscReal :: dstpd_up_ddenup, dstpd_dn_ddendn
  PetscReal :: dsatup, dsatdn
  PetscReal :: delta_X_whatever_dxmolup, delta_X_whatever_dxmoldn
  PetscReal :: dxmass_air_up_dxmol_air_up, dxmass_air_dn_dxmol_air_dn
  PetscReal :: dtot_mole_flux_dstpd, dtot_mole_flux_ddeltaX
  PetscReal :: dtot_mole_flux_ddenave
  PetscReal :: diffusion_scale
  PetscReal :: ddiffusion_coef_dTup, ddiffusion_coef_dTdn
  PetscReal :: ddiffusion_coef_dpup, ddiffusion_coef_dpdn
  PetscReal :: dtot_mole_flux_ddiffusion_coef
  PetscReal :: dstpd_ave_over_dist_dstpd_up, dstpd_ave_over_dist_dstpd_dn
  
  ! Conduction
  PetscReal :: dkeff_up_dsatlup, dkeff_dn_dsatldn
  PetscReal :: dkeff_ave_dkeffup, dkeff_ave_dkeffdn
  PetscReal :: dheat_flux_ddelta_temp, dheat_flux_dkeff_ave
  
  ! DELETE
  
  PetscReal :: Jlup(3,3), Jldn(3,3)
  PetscReal :: Jgup(3,3), Jgdn(3,3)
  PetscReal :: Jcup(3,3), Jcdn(3,3)

  PetscReal :: energy_flux
  PetscReal :: liq_sat, gas_sat, hyd_sat
  PetscReal :: v_sed
  PetscInt  :: gid, lid, hid

  lid = 1
  gid = 2
  hid = 3
 
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
#if 0
!TODO(geh): remove for now
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_up%fracture)) then
    call FracturePermEvaluate(material_auxvar_up,perm_up,temp_perm_up, &
                              dummy_dperm_up,dist)
    perm_up = temp_perm_up
  endif
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif
#endif
  
  if (associated(klinkenberg)) then
    perm_ave_over_dist(1) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
    perm3(:) = perm_up                                        
    perm3 = klinkenberg%Evaluate(perm3, &
                                 hyd_auxvar_up%pres(option%gas_phase))
    temp_perm_up = perm3(1)
    perm3(:) = perm_dn                                        
    perm3 = klinkenberg%Evaluate(perm3, &
                                 hyd_auxvar_dn%pres(option%gas_phase))
    temp_perm_dn = perm3(1)
    perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
                            (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
  else
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  endif
      
  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0  
  
  v_darcy = 0.d0

#ifdef CONVECTION
#ifdef LIQUID_DARCY_FLUX
  iphase = LIQUID_PHASE
  if (hyd_auxvar_up%mobility(iphase) + &
      hyd_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = HydrateAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           hyd_auxvar_up%den_kg, &
                                           hyd_auxvar_dn%den_kg, &
                                           ddensity_kg_ave_dden_kg_up, &
                                           ddensity_kg_ave_dden_kg_dn)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = hyd_auxvar_up%pres(iphase) - &
                     hyd_auxvar_dn%pres(iphase) + &
                     gravity_term
    up_scale = 0.d0
    dn_scale = 0.d0
    upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                             .not.analytical_derivatives, &
                             count_upwind_direction_flip_, &
                             liq_upwind_flip_count_by_res, &
                             liq_upwind_flip_count_by_jac)
    if (upwind) then
      up_scale = 1.d0
      mobility = hyd_auxvar_up%mobility(iphase)
      xmol(:) = hyd_auxvar_up%xmol(:,iphase)
      uH = hyd_auxvar_up%H(iphase)
    else
      dn_scale = 1.d0
      mobility = hyd_auxvar_dn%mobility(iphase)
      xmol(:) = hyd_auxvar_dn%xmol(:,iphase)
      uH = hyd_auxvar_dn%H(iphase)
    endif      

    if (mobility > floweps ) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = HydrateAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          hyd_auxvar_up%den, &
                                          hyd_auxvar_dn%den, &
                                          ddensity_ave_dden_up, &
                                          ddensity_ave_dden_dn)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      tot_mole_flux = q*density_ave
      tot_mole_flux_ddel_pressure = perm_ave_over_dist(iphase) * &
                                       mobility * area * density_ave
      ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
      !                                 xmol[kmol comp/kmol phase]
      wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
      air_mole_flux = tot_mole_flux * xmol(air_comp_id)
      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
    endif                   
  endif
#endif
#ifdef GAS_DARCY_FLUX
  iphase = GAS_PHASE
  if (hyd_auxvar_up%mobility(iphase) + &
      hyd_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = HydrateAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           hyd_auxvar_up%den_kg, &
                                           hyd_auxvar_dn%den_kg, &
                                           ddensity_kg_ave_dden_kg_up, &
                                           ddensity_kg_ave_dden_kg_dn)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = hyd_auxvar_up%pres(iphase) - &
                     hyd_auxvar_dn%pres(iphase) + &
                     gravity_term
    ! if a gas phase does not exist on either side of the connection, the gas
    ! phase properties from the opposite side are used.
    up_scale = 0.d0
    dn_scale = 0.d0
    upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                             .not.analytical_derivatives, &
                             count_upwind_direction_flip_, &
                             gas_upwind_flip_count_by_res, &
                             gas_upwind_flip_count_by_jac)
    if (upwind) then
      up_scale = 1.d0
      mobility = hyd_auxvar_up%mobility(iphase)
      xmol(:) = hyd_auxvar_up%xmol(:,iphase)
      uH = hyd_auxvar_up%H(iphase)
    else
      dn_scale = 1.d0
      mobility = hyd_auxvar_dn%mobility(iphase)
      xmol(:) = hyd_auxvar_dn%xmol(:,iphase)
      uH = hyd_auxvar_dn%H(iphase)
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = HydrateAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          hyd_auxvar_up%den, &
                                          hyd_auxvar_dn%den, &
                                          ddensity_ave_dden_up, &
                                          ddensity_ave_dden_dn)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      tot_mole_flux = q*density_ave
      tot_mole_flux_ddel_pressure = perm_ave_over_dist(iphase) * &
                                       mobility * area * density_ave      
      ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
      !                                 xmol[kmol comp/kmol phase]
      wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
      air_mole_flux = tot_mole_flux * xmol(air_comp_id)
      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
    endif               
  endif
#endif  
  ! CONVECTION
#endif

  ! Sedimentation flux: hydrate
  
  ! q[m^3/sec] = sedimentation velocity[m/sec] * area[m^2]
  ! need to make sure this has a direction, so condition upon gravity?
  ! sedimentation and methanogenesis are linked right now
  if (HYDRATE_WITH_SEDIMENTATION .and. HYDRATE_WITH_METHANOGENESIS) then
    v_sed = methanogenesis%omega
    dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

    if (dabs(dist_gravity) > 0.d0) then
      wat_mole_flux = 0.d0
      air_mole_flux = 0.d0

      q = v_sed * area
  
      upwind = dist_gravity > 0.d0

      if (dist_gravity < 0.d0) q = -q    

      if (upwind) then

        up_scale = 1.d0
      
        hyd_sat = hyd_auxvar_up%sat(hid)
        gas_sat = min(hyd_auxvar_up%sat(gid),hyd_auxvar_up%srg)
        liq_sat = min(hyd_auxvar_up%sat(lid),hyd_auxvar_up%srl)

        wat_mole_flux = hyd_auxvar_up%den(lid)*hyd_auxvar_up%xmol(wat_comp_id, &
                        lid)*liq_sat
        wat_mole_flux = wat_mole_flux + hyd_auxvar_up%den(gid)*hyd_auxvar_up%&
                        xmol(wat_comp_id,gid)*gas_sat
        wat_mole_flux = wat_mole_flux + hyd_auxvar_up%den(hid)*hyd_auxvar_up%&
                        xmol(wat_comp_id,hid)*hyd_sat
        wat_mole_flux = q  * wat_mole_flux

        air_mole_flux = hyd_auxvar_up%den(lid)*hyd_auxvar_up%xmol(air_comp_id, &
                        lid)*liq_sat
        air_mole_flux = air_mole_flux + hyd_auxvar_up%den(gid)*hyd_auxvar_up% &
                        xmol(air_comp_id,gid)*gas_sat
        air_mole_flux = air_mole_flux + hyd_auxvar_up%den(hid)*hyd_auxvar_up% &
                        xmol(air_comp_id,hid)*hyd_sat
        air_mole_flux = q  * air_mole_flux


        energy_flux = q * hyd_auxvar_up%effective_porosity* &
           (hyd_auxvar_up%den(lid) * hyd_auxvar_up%H(lid) * &
           liq_sat + hyd_auxvar_up%den(gid) * hyd_auxvar_up%H(gid) * gas_sat + &
           hyd_auxvar_up%den(hid) * hyd_auxvar_up%H(hid) * hyd_sat) 

      else
        dn_scale = 1.d0

        hyd_sat = hyd_auxvar_dn%sat(hid)
        gas_sat = min(hyd_auxvar_dn%sat(gid),hyd_auxvar_dn%srg)
        liq_sat = min(hyd_auxvar_dn%sat(lid),hyd_auxvar_dn%srl)

        wat_mole_flux = hyd_auxvar_dn%den(lid)*hyd_auxvar_dn%xmol(wat_comp_id, &
                        lid)*liq_sat
        wat_mole_flux = wat_mole_flux + hyd_auxvar_dn%den(gid)*hyd_auxvar_dn%&
                        xmol(wat_comp_id,gid)*gas_sat
        wat_mole_flux = wat_mole_flux + hyd_auxvar_dn%den(hid)*hyd_auxvar_dn%&
                        xmol(wat_comp_id,hid)*hyd_sat
        wat_mole_flux = q  * wat_mole_flux

        air_mole_flux = hyd_auxvar_dn%den(lid)*hyd_auxvar_dn%xmol(air_comp_id, &
                        lid)*liq_sat
        air_mole_flux = air_mole_flux + hyd_auxvar_dn%den(gid)*hyd_auxvar_dn% &
                        xmol(air_comp_id,gid)*gas_sat
        air_mole_flux = air_mole_flux + hyd_auxvar_dn%den(hid)*hyd_auxvar_dn% &
                        xmol(air_comp_id,hid)*hyd_sat
        air_mole_flux = q  * air_mole_flux


        energy_flux = q * hyd_auxvar_dn%effective_porosity * &
                      (hyd_auxvar_dn%den(lid) * hyd_auxvar_dn%H(lid) * &
                      liq_sat + hyd_auxvar_dn%den(gid) * hyd_auxvar_dn%H(gid)*&
                      gas_sat + hyd_auxvar_dn%den(hid) * hyd_auxvar_dn%H(hid)*&
                      hyd_sat)

      endif

      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      Res(energy_id) = Res(energy_id) + energy_flux

    endif
  endif
#ifdef DIFFUSION
  if (.not.hydrate_immiscible) then
  ! add in gas component diffusion in gas and liquid phases
!#if 0
#ifdef LIQUID_DIFFUSION  
  iphase = LIQUID_PHASE
  sat_up = hyd_auxvar_up%sat(iphase)
  sat_dn = hyd_auxvar_dn%sat(iphase)
  dsatup = 1.d0
  dsatdn = 1.d0
  ! by changing #if 1 -> 0, gas component is allowed to diffuse in liquid
  ! phase even if the phase does not exist.
#if 1
  if (sqrt(sat_up*sat_dn) > eps) then
#else
  if (sat_up > eps .or. sat_dn > eps) then
    ! for now, if liquid state neighboring gas, we allow for minute
    ! diffusion in liquid phase.
    if (iphase == option%liquid_phase) then
      if ((sat_up > eps .or. sat_dn > eps)) then
        ! sat_up = max(sat_up,eps)
        if (sat_up < eps) then
          sat_up = eps
          dsatup = 0.d0
        endif
        ! sat_dn = max(sat_dn,eps)
        if (sat_dn < eps) then
          sat_dn = eps
          dsatdn = 0.d0
        endif
      endif
    endif
#endif
    if (hydrate_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_up = hyd_auxvar_up%den(iphase)
      den_dn = hyd_auxvar_dn%den(iphase)
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_up = 1.d0
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = HydrateAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      hyd_auxvar_up%den, &
                                      hyd_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_up = sat_up*material_auxvar_up%tortuosity* &
              hyd_auxvar_up%effective_porosity*den_up
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              hyd_auxvar_dn%effective_porosity*den_dn
              
    dstpd_up_dporup = stpd_up / hyd_auxvar_up%effective_porosity
    dstpd_dn_dpordn = stpd_dn / hyd_auxvar_dn%effective_porosity
    dstpd_up_dsatup = stpd_up / sat_up
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_up_ddenup = tempreal * stpd_up / den_up
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    tempreal = stpd_up*dist_dn+stpd_dn*dist_up
    stpd_ave_over_dist = stpd_up*stpd_dn/tempreal
    dstpd_ave_over_dist_dstpd_up = (stpd_dn-stpd_ave_over_dist*dist_dn)/tempreal
    dstpd_ave_over_dist_dstpd_dn = (stpd_up-stpd_ave_over_dist*dist_up)/tempreal
    
    if (hydrate_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = hyd_auxvar_up%xmol(air_comp_id,iphase) - &
                   hyd_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmolup = 1.d0
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = hyd_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = hyd_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      dxmass_air_up_dxmol_air_up = (fmw_comp(2) - xmass_air_up * (fmw_comp(2) - fmw_comp(1))) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      dxmass_air_dn_dxmol_air_dn = (fmw_comp(2) - xmass_air_dn * (fmw_comp(2) - fmw_comp(1))) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmolup = 1.d0 * dxmass_air_up_dxmol_air_up
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             hydrate_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
    
  endif
#endif
!#if 0
#ifdef GAS_DIFFUSION
  iphase = GAS_PHASE
  sat_up = hyd_auxvar_up%sat(iphase)
  sat_dn = hyd_auxvar_dn%sat(iphase)
  !geh: i am not sure why both of these conditionals were included.  seems
  !     like the latter would never be false.
  if (sqrt(sat_up*sat_dn) > eps) then
    dsatup = 1.d0
    dsatdn = 1.d0
    if (hydrate_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_up = hyd_auxvar_up%den(iphase)
      den_dn = hyd_auxvar_dn%den(iphase)
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_up = 1.d0
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = HydrateAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      hyd_auxvar_up%den, &
                                      hyd_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_up = sat_up*material_auxvar_up%tortuosity* &
              hyd_auxvar_up%effective_porosity*den_up
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              hyd_auxvar_dn%effective_porosity*den_dn
              
    dstpd_up_dporup = stpd_up / hyd_auxvar_up%effective_porosity
    dstpd_dn_dpordn = stpd_dn / hyd_auxvar_dn%effective_porosity
    dstpd_up_dsatup = stpd_up / sat_up
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_up_ddenup = tempreal * stpd_up / den_up
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    tempreal = stpd_up*dist_dn+stpd_dn*dist_up
    stpd_ave_over_dist = stpd_up*stpd_dn/tempreal
    dstpd_ave_over_dist_dstpd_up = (stpd_dn-stpd_ave_over_dist*dist_dn)/tempreal
    dstpd_ave_over_dist_dstpd_dn = (stpd_up-stpd_ave_over_dist*dist_up)/tempreal
    
    if (hydrate_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = hyd_auxvar_up%xmol(air_comp_id,iphase) - &
                   hyd_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmolup = 1.d0
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = hyd_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = hyd_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      dxmass_air_up_dxmol_air_up = (fmw_comp(2) - xmass_air_up * (fmw_comp(2) - fmw_comp(1))) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      dxmass_air_dn_dxmol_air_dn = (fmw_comp(2) - xmass_air_dn * (fmw_comp(2) - fmw_comp(1))) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmolup = 1.d0 * dxmass_air_up_dxmol_air_up
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    ! need to account for multiple phases
    ! Eq. 1.9b.  The gas density is added below
    if (hydrate_temp_dep_gas_air_diff) then
      temp_ave = 0.5d0*(hyd_auxvar_up%temp+hyd_auxvar_dn%temp)
      pressure_ave = 0.5d0*(hyd_auxvar_up%pres(iphase)+ &
                            hyd_auxvar_dn%pres(iphase))
      tempreal = (temp_ave+273.15d0)/273.15d0
      diffusion_scale = tempreal**1.8d0 * 101325.d0 / pressure_ave
                             ! 0.9d0 = 0.5 * 1.8
      ddiffusion_coef_dTup = 0.9d0 * diffusion_scale / (tempreal * 273.15d0)
      ddiffusion_coef_dTdn = ddiffusion_coef_dTup
      ddiffusion_coef_dpup = -1.d0 * diffusion_scale / pressure_ave * 0.5d0
      ddiffusion_coef_dpdn = ddiffusion_coef_dpup
    else
      diffusion_scale = 1.d0
      ddiffusion_coef_dTup = 0.d0
      ddiffusion_coef_dTdn = 0.d0
      ddiffusion_coef_dpup = 0.d0
      ddiffusion_coef_dpdn = 0.d0
    endif
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             diffusion_scale * &
                             hydrate_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddiffusion_coef = tot_mole_flux / diffusion_scale
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave    
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
  endif
#endif
! DIFFUSION
  endif ! if (.not.hydrate_immiscible)
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_wet-k_dry)
  ! 1 = dry
  ! 2 = wet
  sat_up = hyd_auxvar_up%sat(option%liquid_phase)
  sat_dn = hyd_auxvar_dn%sat(option%liquid_phase)
  
  call HydrateCompositeThermalCond(material_auxvar_up%porosity, &
               hyd_auxvar_up%sat,thermal_conductivity_up(1), &
               thermal_conductivity_dn(2),k_eff_up)
  call HydrateCompositeThermalCond(material_auxvar_dn%porosity, &
              hyd_auxvar_dn%sat,thermal_conductivity_dn(1), &
              thermal_conductivity_dn(2), k_eff_dn)
  
  if (k_eff_up > 0.d0 .or. k_eff_dn > 0.d0) then
    tempreal = k_eff_up*dist_dn+k_eff_dn*dist_up
    k_eff_ave = k_eff_up*k_eff_dn/tempreal
    dkeff_ave_dkeffup = (k_eff_dn-k_eff_ave*dist_dn)/tempreal
    dkeff_ave_dkeffdn = (k_eff_up-k_eff_ave*dist_up)/tempreal
  else
    k_eff_ave = 0.d0
    dkeff_ave_dkeffup = 0.d0
    dkeff_ave_dkeffdn = 0.d0
  endif
  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = hyd_auxvar_up%temp - hyd_auxvar_dn%temp
  dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
  heat_flux = dheat_flux_ddelta_temp * delta_temp
  dheat_flux_dkeff_ave = area * 1.d-6 * delta_temp
  ! MJ/s or MW
  Res(energy_id) = Res(energy_id) + heat_flux
  
! CONDUCTION
#endif

end subroutine HydrateFlux

! ************************************************************************** !

subroutine HydrateBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         hyd_auxvar_up,global_auxvar_up, &
                         hyd_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         thermal_conductivity_dn, &
                         area,dist,upwind_direction_, &
                         methanogenesis, &
                         hydrate_parameter, &
                         option,v_darcy,Res,J, &
                         analytical_derivatives, &
                         update_upwind_direction_, &
                         count_upwind_direction_flip_, &
                         debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Option_module                              
  use Material_Aux_class
  use Fracture_module
  use Klinkenberg_module
  use Upwind_Direction_module

  implicit none
  
  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(HYDRATE_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(hydrate_auxvar_type) :: hyd_auxvar_up, hyd_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(methanogenesis_type) :: methanogenesis
  type(hydrate_parameter_type) :: hydrate_parameter
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: J(3,3)
  PetscBool :: analytical_derivatives
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_
  PetscBool :: debug_connection
  
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: xmol(option%nflowspec)  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_xmol, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, q 
  PetscReal :: tot_mole_flux
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_dn, k_eff_ave, heat_flux
  PetscReal :: boundary_pressure
  PetscReal :: xmass_air_up, xmass_air_dn, delta_xmass  
  PetscReal :: xmol_air_up, xmol_air_dn
  PetscReal :: tempreal
  PetscReal :: delta_X_whatever
  PetscReal :: wat_mole_flux, air_mole_flux
  PetscBool :: upwind

  ! Darcy flux
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dpaup, ddelta_pressure_dpadn
  PetscReal :: ddelta_pressure_dTup, ddelta_pressure_dTdn
  PetscReal :: dv_darcy_ddelta_pressure
  PetscReal :: dv_darcy_dmobility
  
  PetscReal :: up_scale, dn_scale
  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dtot_mole_flux_dp, dtot_mole_flux_dT, dtot_mole_flux_dsatg
  PetscReal :: dpl_dsatg
  PetscReal :: ddelta_pressure_pl
  PetscReal :: tot_mole_flux_ddel_pressure, tot_mole_flux_dmobility
  PetscReal :: xmol_bool

  ! Diffusion
  PetscReal :: stpd_dn
  PetscReal :: dstpd_dn_dpordn
  PetscReal :: dstpd_dn_dsatdn
  PetscReal :: dstpd_dn_ddendn
  PetscReal :: dsatdn
  PetscReal :: delta_X_whatever_dxmoldn
  PetscReal :: dxmass_air_dn_dxmol_air_dn
  PetscReal :: dtot_mole_flux_dstpd, dtot_mole_flux_ddeltaX
  PetscReal :: dtot_mole_flux_ddenave
  PetscReal :: diffusion_scale
  PetscReal :: ddiffusion_coef_dTdn
  PetscReal :: ddiffusion_coef_dpdn
  PetscReal :: dtot_mole_flux_ddiffusion_coef
  PetscReal :: dstpd_ave_over_dist_dstpd_dn
  PetscReal :: pressure_ave
  PetscReal :: perm3(3)
  
  ! Conduction
  PetscReal :: dkeff_up_dsatlup, dkeff_dn_dsatldn
  PetscReal :: dkeff_ave_dkeffup, dkeff_ave_dkeffdn
  PetscReal :: dheat_flux_ddelta_temp, dheat_flux_dkeff_ave
  
  ! DELETE
  
  PetscReal :: Jl(3,3)
  PetscReal :: Jg(3,3)
  PetscReal :: Jc(3,3)
  
  PetscInt :: idof
  
  PetscReal :: temp_perm_dn
  PetscReal :: dummy_dperm_dn

  PetscReal :: energy_flux
  PetscReal :: liq_sat, gas_sat, hyd_sat
  PetscReal :: v_sed
  PetscInt  :: gid, lid, hid

  lid = 1
  gid = 2
  hid = 3
 
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id

  Res = 0.d0
  J = 0.d0
  v_darcy = 0.d0  

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

#if 0
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_dn%fracture)) then
    call FracturePermEvaluate(material_auxvar_dn,perm_dn,temp_perm_dn, &
                              dummy_dperm_dn,dist)
    perm_dn = temp_perm_dn
  endif  
#endif
  
  if (associated(klinkenberg)) then
    perm_dn_adj(1) = perm_dn
    perm3(:) = perm_dn                                        
    perm3 = klinkenberg%Evaluate(perm3, &
                                 hyd_auxvar_dn%pres(option%gas_phase))
    perm_dn_adj(2) = perm3(1)
  else
    perm_dn_adj(:) = perm_dn
  endif
  
#ifdef CONVECTION  
#ifdef LIQUID_DARCY_FLUX
  iphase = LIQUID_PHASE
  mobility = 0.d0
  xmol_bool = 1.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      if (hyd_auxvar_up%mobility(iphase) + &
          hyd_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(HYDRATE_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(HYDRATE_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        boundary_pressure = hyd_auxvar_up%pres(iphase)
        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == G_STATE) then
          ! the idea here is to accommodate a free surface boundary
          ! face.  this will not work for an interior grid cell as
          ! there should be capillary pressure in force.
          boundary_pressure = hyd_auxvar_up%pres(option%gas_phase)
        endif
        density_kg_ave = HydrateAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                hyd_auxvar_up%den_kg, &
                                                hyd_auxvar_dn%den_kg, &
                                                ddensity_kg_ave_dden_kg_up, &
                                                ddensity_kg_ave_dden_kg_dn)
        ddensity_kg_ave_dden_kg_up = 0.d0 ! always
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          hyd_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              hyd_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
          endif
        endif
        dn_scale = 0.d0
        upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                                 .not.analytical_derivatives, &
                                 count_upwind_direction_flip_, &
                                 liq_bc_upwind_flip_count_by_res, &
                                 liq_bc_upwind_flip_count_by_jac)
        if (upwind) then
          mobility = hyd_auxvar_up%mobility(iphase)
          xmol(:) = hyd_auxvar_up%xmol(:,iphase)
          uH = hyd_auxvar_up%H(iphase)
        else
          dn_scale = 1.d0        
          mobility = hyd_auxvar_dn%mobility(iphase)
          xmol(:) = hyd_auxvar_dn%xmol(:,iphase)
          uH = hyd_auxvar_dn%H(iphase)
        endif      
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = HydrateAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            hyd_auxvar_up%den, &
                                            hyd_auxvar_dn%den, &
                                            ddensity_ave_dden_up, &
                                            ddensity_ave_dden_dn)    
        ddensity_ave_dden_up = 0.d0 ! always
        dv_darcy_dmobility = perm_ave_over_dist * delta_pressure
      endif
    case(NEUMANN_BC)
      xmol_bool = 0.d0
      dv_darcy_ddelta_pressure = 0.d0
      dv_darcy_dmobility = 0.d0
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      ddelta_pressure_dpdn = 0.d0
      ddelta_pressure_dTdn = 0.d0
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(HYDRATE_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(HYDRATE_GAS_FLUX_INDEX)
      end select
      xmol = 0.d0
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      xmol(iphase) = 1.d0
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = hyd_auxvar_up%den(iphase)
          uH = hyd_auxvar_up%H(iphase)
        else 
          dn_scale = 1.d0
          density_ave = hyd_auxvar_dn%den(iphase)
          uH = hyd_auxvar_dn%H(iphase)
          ddensity_ave_dden_dn = 1.d0
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in HydrateBCFlux phase loop.'
      call PrintErrMsg(option)
  end select
  if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    tot_mole_flux_ddel_pressure = dv_darcy_ddelta_pressure * area * &
                                  density_ave
    tot_mole_flux_dmobility = dv_darcy_dmobility * area * density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
    air_mole_flux = tot_mole_flux * xmol(air_comp_id)
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
    Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
  endif                   
#endif
#ifdef GAS_DARCY_FLUX
  iphase = GAS_PHASE
  mobility = 0.d0
  xmol_bool = 1.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      if (hyd_auxvar_up%mobility(iphase) + &
          hyd_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(HYDRATE_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(HYDRATE_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        boundary_pressure = hyd_auxvar_up%pres(iphase)
        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == G_STATE) then
          ! the idea here is to accommodate a free surface boundary
          ! face.  this will not work for an interior grid cell as
          ! there should be capillary pressure in force.
          boundary_pressure = hyd_auxvar_up%pres(option%gas_phase)
        endif
        density_kg_ave = HydrateAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                hyd_auxvar_up%den_kg, &
                                                hyd_auxvar_dn%den_kg, &
                                                ddensity_kg_ave_dden_kg_up, &
                                                ddensity_kg_ave_dden_kg_dn)
        ddensity_kg_ave_dden_kg_up = 0.d0 ! always
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          hyd_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              hyd_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
          endif
        endif
        dn_scale = 0.d0
        ! don't expect the derivative to match precisely at delta_pressure = 0
        ! due to potential switch in direction for numerically perturbed
        ! residual
        upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                                 .not.analytical_derivatives, &
                                 count_upwind_direction_flip_, &
                                 gas_bc_upwind_flip_count_by_res, &
                                 gas_bc_upwind_flip_count_by_jac)
        if (upwind) then
          mobility = hyd_auxvar_up%mobility(iphase)
          xmol(:) = hyd_auxvar_up%xmol(:,iphase)
          uH = hyd_auxvar_up%H(iphase)
        else
          dn_scale = 1.d0        
          mobility = hyd_auxvar_dn%mobility(iphase)
          xmol(:) = hyd_auxvar_dn%xmol(:,iphase)
          uH = hyd_auxvar_dn%H(iphase)
        endif  
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = HydrateAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            hyd_auxvar_up%den, &
                                            hyd_auxvar_dn%den, &
                                            ddensity_ave_dden_up, &
                                            ddensity_ave_dden_dn)    
        ddensity_ave_dden_up = 0.d0 ! always
        dv_darcy_dmobility = perm_ave_over_dist * delta_pressure
      endif
    case(NEUMANN_BC)
      xmol_bool = 0.d0
      dv_darcy_ddelta_pressure = 0.d0
      dv_darcy_dmobility = 0.d0
      ddensity_ave_dden_up = 0.d0 ! always
      ddensity_ave_dden_dn = 0.d0
      ddelta_pressure_dpdn = 0.d0
      ddelta_pressure_dpadn = 0.d0
      ddelta_pressure_dTdn = 0.d0      
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(HYDRATE_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(HYDRATE_GAS_FLUX_INDEX)
      end select
      xmol = 0.d0
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      xmol(iphase) = 1.d0
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = hyd_auxvar_up%den(iphase)
          uH = hyd_auxvar_up%H(iphase)
        else 
          dn_scale = 1.d0
          density_ave = hyd_auxvar_dn%den(iphase)
          uH = hyd_auxvar_dn%H(iphase)
          ddensity_ave_dden_dn = 1.d0
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in HydrateBCFlux phase loop.'
      call PrintErrMsg(option)
  end select

  if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    tot_mole_flux_ddel_pressure = dv_darcy_ddelta_pressure * area * &
                                  density_ave
    tot_mole_flux_dmobility = dv_darcy_dmobility * area * density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    wat_mole_flux = tot_mole_flux * xmol(wat_comp_id)
    air_mole_flux = tot_mole_flux * xmol(air_comp_id)
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
    Res(energy_id) = Res(energy_id) + tot_mole_flux * uH
  endif                   
#endif  
! CONVECTION
#endif

  ! Sedimentation flux: hydrate

  ! q[m^3/sec] = sedimentation velocity[m/sec] * area[m^2]
  ! need to make sure this has a direction, so condition upon gravity
  ! sedimentation and methanogenesis are linked right now.
  if (HYDRATE_WITH_SEDIMENTATION .and. HYDRATE_WITH_METHANOGENESIS) then
    v_sed = methanogenesis%omega
    dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

    if (dabs(dist_gravity) > 0.d0) then
      wat_mole_flux = 0.d0
      air_mole_flux = 0.d0

      q = v_sed * area

      upwind = dist_gravity > 0.d0

      if (dist_gravity < 0.d0) q = -q

      hyd_sat = hyd_auxvar_dn%sat(hid)
      gas_sat = min(hyd_auxvar_dn%sat(gid),hyd_auxvar_dn%srg)
      liq_sat = min(hyd_auxvar_dn%sat(lid),hyd_auxvar_dn%srl)

      wat_mole_flux = hyd_auxvar_dn%den(lid)*hyd_auxvar_dn%xmol(wat_comp_id, &
                        lid)*liq_sat
      wat_mole_flux = wat_mole_flux + hyd_auxvar_dn%den(gid)*hyd_auxvar_dn%&
                        xmol(wat_comp_id,gid)*gas_sat
      wat_mole_flux = wat_mole_flux + hyd_auxvar_dn%den(hid)*hyd_auxvar_dn%&
                        xmol(wat_comp_id,hid)*hyd_sat
      wat_mole_flux = q  * wat_mole_flux

      air_mole_flux = hyd_auxvar_dn%den(lid)*hyd_auxvar_dn%xmol(air_comp_id, &
                      lid)*liq_sat
      air_mole_flux = air_mole_flux + hyd_auxvar_dn%den(gid)*hyd_auxvar_dn% &
                      xmol(air_comp_id,gid)*gas_sat
      air_mole_flux = air_mole_flux + hyd_auxvar_dn%den(hid)*hyd_auxvar_dn% &
                      xmol(air_comp_id,hid)*hyd_sat
      air_mole_flux = q  * air_mole_flux


    ! MAN: need to mult by phi?
      energy_flux = q*hyd_auxvar_dn%effective_porosity*(hyd_auxvar_dn%den(lid) * &
                     hyd_auxvar_dn%H(lid) * liq_sat + &
                     hyd_auxvar_dn%den(gid) * hyd_auxvar_dn%H(gid) * gas_sat + &
                     hyd_auxvar_dn%den(hid) * hyd_auxvar_dn%H(hid) * hyd_sat)

      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      Res(energy_id) = Res(energy_id) + energy_flux

    endif
  endif

  
#ifdef DIFFUSION
  if (.not.hydrate_immiscible) then
#ifdef LIQUID_DIFFUSION  
  iphase = LIQUID_PHASE
  dsatdn = 1.d0
  ! diffusion all depends upon the downwind cell.  phase diffusion only
  ! occurs if a phase exists in both auxvars (boundary and internal) or
  ! a liquid phase exists in the internal cell. so, one could say that
  ! liquid diffusion always exists as the internal cell has a liquid phase,
  ! but gas phase diffusion only occurs if the internal cell has a gas
  ! phase.
  sat_dn = hyd_auxvar_dn%sat(iphase)
  if (sat_dn > eps .and. ibndtype(iphase) /= NEUMANN_BC) then
    if (hydrate_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_dn = hyd_auxvar_dn%den(iphase)
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      density_ave = HydrateAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      hyd_auxvar_up%den, &
                                      hyd_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ddensity_ave_dden_up = 0.d0
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              hyd_auxvar_dn%effective_porosity*den_dn
              
    dstpd_dn_dpordn = stpd_dn / hyd_auxvar_dn%effective_porosity
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    dstpd_ave_over_dist_dstpd_dn = 1.d0 / dist(0)
    stpd_ave_over_dist = stpd_dn * dstpd_ave_over_dist_dstpd_dn
    
    if (hydrate_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = hyd_auxvar_up%xmol(air_comp_id,iphase) - &
                   hyd_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = hyd_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = hyd_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      dxmass_air_dn_dxmol_air_dn = (fmw_comp(2) - xmass_air_dn * (fmw_comp(2) - fmw_comp(1))) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             hydrate_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
  endif
#endif
#ifdef GAS_DIFFUSION
  iphase = GAS_PHASE
  sat_dn = hyd_auxvar_dn%sat(iphase)
  !geh: i am not sure why both of these conditionals were included.  seems
  !     like the latter would never be false.
  if (sat_dn > eps .and. ibndtype(iphase) /= NEUMANN_BC) then
    dsatdn = 1.d0
    if (hydrate_harmonic_diff_density) then
      ! density_ave in this case is not used.
      density_ave = 1.d0
      den_dn = hyd_auxvar_dn%den(iphase)
      ddensity_ave_dden_dn = 0.d0
      tempreal = 1.d0
    else
      ! den_up and den_dn are not used in this case
      den_dn = 1.d0
      ! we use upstream weighting when iphase is not equal, otherwise
      ! arithmetic with 50/50 weighting
      !TODO(geh): why are we averaging density here?
      density_ave = HydrateAverageDensity(iphase, &
                                      global_auxvar_up%istate, &
                                      global_auxvar_dn%istate, &
                                      hyd_auxvar_up%den, &
                                      hyd_auxvar_dn%den, &
                                      ddensity_ave_dden_up, &
                                      ddensity_ave_dden_dn)
      ddensity_ave_dden_up = 0.d0                                      
      ! used to zero out derivative below
      tempreal = 0.d0
    endif
    stpd_dn = sat_dn*material_auxvar_dn%tortuosity* &
              hyd_auxvar_dn%effective_porosity*den_dn
              
    dstpd_dn_dpordn = stpd_dn / hyd_auxvar_dn%effective_porosity
    dstpd_dn_dsatdn = stpd_dn / sat_dn
    dstpd_dn_ddendn = tempreal * stpd_dn / den_dn
    ! units = [mole/m^4 bulk]
    dstpd_ave_over_dist_dstpd_dn = 1.d0 / dist(0)
    stpd_ave_over_dist = stpd_dn * dstpd_ave_over_dist_dstpd_dn    
    
    if (hydrate_diffuse_xmol) then ! delta of mole fraction
      delta_xmol = hyd_auxvar_up%xmol(air_comp_id,iphase) - &
                   hyd_auxvar_dn%xmol(air_comp_id,iphase)
      delta_X_whatever = delta_xmol
      delta_X_whatever_dxmoldn = -1.d0
    else ! delta of mass fraction
      xmol_air_up = hyd_auxvar_up%xmol(air_comp_id,iphase)
      xmol_air_dn = hyd_auxvar_dn%xmol(air_comp_id,iphase)
      tempreal = (xmol_air_up*fmw_comp(2) + (1.d0-xmol_air_up)*fmw_comp(1))
      xmass_air_up = xmol_air_up*fmw_comp(2) / tempreal
      tempreal = (xmol_air_dn*fmw_comp(2) + (1.d0-xmol_air_dn)*fmw_comp(1))
      xmass_air_dn = xmol_air_dn*fmw_comp(2) / tempreal
      delta_xmass = xmass_air_up - xmass_air_dn
      delta_X_whatever = delta_xmass
      delta_X_whatever_dxmoldn = -1.d0 * dxmass_air_dn_dxmol_air_dn
    endif
    ! need to account for multiple phases
    ! Eq. 1.9b.  The gas density is added below
    if (hydrate_temp_dep_gas_air_diff) then
      temp_ave = 0.5d0*(hyd_auxvar_up%temp+hyd_auxvar_dn%temp)
      pressure_ave = 0.5d0*(hyd_auxvar_up%pres(iphase)+ &
                            hyd_auxvar_dn%pres(iphase))
      tempreal = (temp_ave+273.15d0)/273.15d0
      diffusion_scale = tempreal**1.8d0 * 101325.d0 / pressure_ave
                             ! 0.9d0 = 0.5 * 1.8
      ddiffusion_coef_dTdn = 0.9d0 * diffusion_scale / (tempreal * 273.15d0)
      ddiffusion_coef_dpdn = -1.d0 * diffusion_scale / pressure_ave * 0.5d0
    else
      diffusion_scale = 1.d0
      ddiffusion_coef_dTdn = 0.d0
      ddiffusion_coef_dpdn = 0.d0
    endif
    ! units = mole/sec
    dtot_mole_flux_ddeltaX = density_ave * stpd_ave_over_dist * &
                             diffusion_scale * &
                             hydrate_parameter%diffusion_coefficient(iphase) * &
                             area
    tot_mole_flux = dtot_mole_flux_ddeltaX * delta_X_whatever
    dtot_mole_flux_dstpd = tot_mole_flux / stpd_ave_over_dist
    dtot_mole_flux_ddiffusion_coef = tot_mole_flux / diffusion_scale
    dtot_mole_flux_ddenave = tot_mole_flux / density_ave    
    Res(wat_comp_id) = Res(wat_comp_id) - tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + tot_mole_flux
  endif
#endif
! DIFFUSION
  endif ! if (.not.hydrate_immiscible)
#endif

#ifdef CONDUCTION
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_wet-k_dry)
  ! 1 = dry
  ! 2 = wet
  heat_flux = 0.d0
  select case (ibndtype(HYDRATE_ENERGY_EQUATION_INDEX))
    case (DIRICHLET_BC)
      sat_dn = hyd_auxvar_dn%sat(option%liquid_phase)
      call HydrateCompositeThermalCond(material_auxvar_dn%porosity, &
                  hyd_auxvar_dn%sat,thermal_conductivity_dn(1), &
                  thermal_conductivity_dn(2),k_eff_dn)

      dkeff_ave_dkeffdn = 1.d0 / dist(0)
      k_eff_ave = k_eff_dn * dkeff_ave_dkeffdn
      ! units:
      ! k_eff = W/K-m = J/s/K-m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = k_eff * delta_temp * area = J/s
      delta_temp = hyd_auxvar_up%temp - hyd_auxvar_dn%temp
      dheat_flux_ddelta_temp = k_eff_ave * area * 1.d-6 ! J/s -> MJ/s
      heat_flux = dheat_flux_ddelta_temp * delta_temp
      dheat_flux_dkeff_ave = area * 1.d-6 * delta_temp
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = auxvars(auxvar_mapping(HYDRATE_ENERGY_FLUX_INDEX)) * area
      dheat_flux_ddelta_temp = 0.d0
      dkeff_dn_dsatldn = 0.d0
      dkeff_ave_dkeffdn = 0.d0
      dheat_flux_dkeff_ave = 0.d0
    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'HydrateBCFlux heat conduction loop.'
      call PrintErrMsg(option)
  end select
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
  
! CONDUCTION
#endif

end subroutine HydrateBCFlux

! ************************************************************************** !

subroutine HydrateSrcSink(option,qsrc,flow_src_sink_type,hyd_auxvar_ss, &
                          hyd_auxvar,global_auxvar,ss_flow_vol_flux, &
                          scale,Res,J,analytical_derivatives,debug_cell)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  type(hydrate_auxvar_type) :: hyd_auxvar,hyd_auxvar_ss
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscReal :: J(option%nflowdof,option%nflowdof)  
  PetscBool :: analytical_derivatives  
  PetscBool :: debug_cell
  
  PetscReal :: qsrc(3)
  PetscInt :: flow_src_sink_type
  PetscReal :: qsrc_mol
  PetscReal :: enthalpy, internal_energy
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscReal :: Jl(option%nflowdof,option%nflowdof)  
  PetscReal :: Jg(option%nflowdof,option%nflowdof)  
  PetscReal :: Je(option%nflowdof,option%nflowdof)  
  PetscReal :: dden_bool
  PetscReal :: hw_dp, hw_dT, ha_dp, ha_dT
  PetscErrorCode :: ierr
  PetscReal :: mob_tot
  PetscInt, parameter :: lid = 1
  PetscInt, parameter :: gid = 2

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  energy_id = option%energy_id
  
  Res = 0.d0
  J = 0.d0
 
  qsrc_mol = 0.d0 
  if (flow_src_sink_type == TOTAL_MASS_RATE_SS) then
    !MAN: this has only been tested for an extraction well. Scales the mass of 
    !water and gas extracted by the mobility ratio.
    mob_tot = hyd_auxvar%mobility(lid) + hyd_auxvar%mobility(gid)
      if (hyd_auxvar%sat(gid) <= 0.d0) then
        ! Water component, liquid phase
        ! kg/s phase to kmol/sec phase
        qsrc_mol = qsrc(wat_comp_id) * hyd_auxvar%den(lid) / hyd_auxvar%den_kg(lid)
        ! kmol/sec phase to kmol/sec component
        qsrc_mol = qsrc_mol * hyd_auxvar%xmol(lid,lid)
      elseif (hyd_auxvar%sat(lid) <= 0.d0) then
        ! Water component, gas phase
        ! kg/s phase to kmol/sec phase
        qsrc_mol = qsrc(wat_comp_id) * hyd_auxvar%den(gid) / hyd_auxvar%den_kg(gid)
        ! kmol/sec phase to kmol/sec component
        qsrc_mol = qsrc_mol * hyd_auxvar%xmol(lid,gid)
      else
        ! Water component, liquid phase
        qsrc_mol = qsrc(wat_comp_id) * hyd_auxvar%den(lid) / &
                   hyd_auxvar%den_kg(lid)*hyd_auxvar%mobility(lid)/mob_tot * &
                   hyd_auxvar%xmol(lid,lid)
        ! Water component, gas phase
        qsrc_mol = qsrc_mol + qsrc(wat_comp_id) * hyd_auxvar%den(gid) / &
                   hyd_auxvar%den_kg(gid)*hyd_auxvar%mobility(gid)/mob_tot * &
                   hyd_auxvar%xmol(lid,gid)
      endif
      
      ss_flow_vol_flux(wat_comp_id) = qsrc_mol/hyd_auxvar%den(wat_comp_id)
      Res(wat_comp_id) = qsrc_mol
      
      if (hyd_auxvar%sat(gid) <= 0.d0) then
        ! Air component, liquid phase
        ! kg/s phase to kmol/sec phase
        qsrc_mol = qsrc(wat_comp_id) * hyd_auxvar%den(lid) / hyd_auxvar%den_kg(lid)
        ! kmol/sec phase to kmol/sec component
        qsrc_mol = qsrc_mol * hyd_auxvar%xmol(gid,lid)
      elseif (hyd_auxvar%sat(lid) <= 0.d0) then
        ! Air component, gas phase
        ! kg/s phase to kmol/sec phase
        qsrc_mol = qsrc(wat_comp_id) * hyd_auxvar%den(gid) / hyd_auxvar%den_kg(gid)
        ! kmol/sec phase to kmol/sec component
        qsrc_mol = qsrc_mol * hyd_auxvar%xmol(gid,gid)
      else
      ! Air component, liquid phase
        qsrc_mol = qsrc(wat_comp_id) * hyd_auxvar%den(lid) / &
                   hyd_auxvar%den_kg(lid)*hyd_auxvar%mobility(lid)/mob_tot * &
                   hyd_auxvar%xmol(gid,lid)
      ! Air component, gas phase
        qsrc_mol = qsrc_mol + qsrc(wat_comp_id) * hyd_auxvar%den(gid) / &
                   hyd_auxvar%den_kg(gid)*hyd_auxvar%mobility(gid)/mob_tot * &
                   hyd_auxvar%xmol(gid,gid)
      endif
   
      ss_flow_vol_flux(air_comp_id) = qsrc_mol/hyd_auxvar%den(air_comp_id)
      Res(air_comp_id) = qsrc_mol

      if (hyd_auxvar%sat(gid) <= 0.d0) then
        Res(energy_id) = qsrc(wat_comp_id) * hyd_auxvar%den(lid) / &
                         hyd_auxvar%den_kg(lid) * hyd_auxvar_ss%h(lid)
      elseif (hyd_auxvar%sat(lid) <= 0.d0) then
        Res(energy_id) = qsrc(wat_comp_id) * hyd_auxvar%den(gid) / &
                         hyd_auxvar%den_kg(gid) * hyd_auxvar_ss%h(gid)
      else
        Res(energy_id) = qsrc(wat_comp_id) * hyd_auxvar%mobility(lid)/mob_tot*&
                         hyd_auxvar%den(lid) / hyd_auxvar%den_kg(lid) * &
                         hyd_auxvar_ss%h(lid)
        Res(energy_id) = Res(energy_id) + qsrc(wat_comp_id) * hyd_auxvar% &
                         mobility(gid)/mob_tot*hyd_auxvar%den(gid) / &
                         hyd_auxvar%den_kg(gid) * hyd_auxvar_ss%h(gid)
      endif
  else
  
#ifdef WATER_SRCSINK
  qsrc_mol = 0.d0
  dden_bool = 0.d0
  select case(flow_src_sink_type)
    case(MASS_RATE_SS)
      qsrc_mol = qsrc(wat_comp_id)/fmw_comp(wat_comp_id) ! kg/sec -> kmol/sec
    case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
      qsrc_mol = qsrc(wat_comp_id)/fmw_comp(wat_comp_id)*scale 
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_mol = qsrc(wat_comp_id)*hyd_auxvar%den(wat_comp_id) ! den = kmol/m^3
      dden_bool = 1.d0
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_mol = qsrc(wat_comp_id)*hyd_auxvar%den(wat_comp_id)*scale
      dden_bool = 1.d0
  end select
  ss_flow_vol_flux(wat_comp_id) = qsrc_mol/hyd_auxvar%den(wat_comp_id)
  Res(wat_comp_id) = qsrc_mol
#endif

#ifdef AIR_SRCSINK
  qsrc_mol = 0.d0
  dden_bool = 0.d0
  select case(flow_src_sink_type)
    case(MASS_RATE_SS)
      qsrc_mol = qsrc(air_comp_id)/fmw_comp(air_comp_id) ! kg/sec -> kmol/sec
    case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
      qsrc_mol = qsrc(air_comp_id)/fmw_comp(air_comp_id)*scale 
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_mol = qsrc(air_comp_id)*hyd_auxvar%den(air_comp_id) ! den = kmol/m^3
      dden_bool = 1.d0
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_mol = qsrc(air_comp_id)*hyd_auxvar%den(air_comp_id)*scale
      dden_bool = 1.d0
  end select
  ss_flow_vol_flux(air_comp_id) = qsrc_mol/hyd_auxvar%den(air_comp_id)
  Res(air_comp_id) = qsrc_mol
#endif
  endif
  
  if (dabs(qsrc(air_comp_id)) < 1.d-40 .and. flow_src_sink_type /= &
          TOTAL_MASS_RATE_SS .and. qsrc(wat_comp_id) < 0.d0) then ! extraction only
    ! Res(1) holds qsrc_mol for water.  If the src/sink value for air is zero,
    ! remove/add the equivalent mole fraction of air in the liquid phase.
    qsrc_mol = Res(wat_comp_id)*hyd_auxvar%xmol(air_comp_id,wat_comp_id)
    Res(air_comp_id) = qsrc_mol
    ss_flow_vol_flux(air_comp_id) = qsrc_mol/hyd_auxvar%den(air_comp_id)
  endif
  
  ! energy units: MJ/sec
  if (size(qsrc) == THREE_INTEGER) then
    if (flow_src_sink_type /= TOTAL_MASS_RATE_SS) then
      if (dabs(qsrc(wat_comp_id)) > 1.d-40) then
        if (associated(hyd_auxvar%d)) then
          hw_dp = hyd_auxvar_ss%d%Hl_pl          
          hw_dT = hyd_auxvar_ss%d%Hl_T
        endif
        enthalpy = hyd_auxvar_ss%h(wat_comp_id)
        ! enthalpy units: MJ/kmol                       ! water component mass
        Res(energy_id) = Res(energy_id) + Res(wat_comp_id) * enthalpy       
        J = J + Je
      endif
      
      if (dabs(qsrc(air_comp_id)) > 1.d-40) then
        ! this is pure air, we use the enthalpy of air, NOT the air/water
        ! mixture in gas
        ! air enthalpy is only a function of temperature
        if (associated(hyd_auxvar%d)) then
          ha_dp = hyd_auxvar_ss%d%Ha_pg
          ha_dT = hyd_auxvar_ss%d%Ha_T
        endif
        
        internal_energy = hyd_auxvar_ss%u(air_comp_id)
        enthalpy = hyd_auxvar_ss%h(air_comp_id)                                 
        ! enthalpy units: MJ/kmol                       ! air component mass
        Res(energy_id) = Res(energy_id) + Res(air_comp_id) * enthalpy
        J = J + Je
      endif
    endif
    Res(energy_id) = Res(energy_id) + qsrc(energy_id)*scale ! MJ/s
    ! no derivative
  endif
  
end subroutine HydrateSrcSink

! ************************************************************************** !

subroutine HydrateAccumDerivative(hyd_auxvar,global_auxvar,material_auxvar, &
                                  z,offset,meth,soil_heat_capacity,option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(hydrate_auxvar_type) :: hyd_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: z,offset
  type(methanogenesis_type), pointer :: meth
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof) 
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow
 
  if (.not. hydrate_central_diff_jacobian) then
    call HydrateAccumulation(hyd_auxvar(ZERO_INTEGER),global_auxvar, &
                           material_auxvar,z,offset,meth,soil_heat_capacity, &
                           option,res,jac,hydrate_analytical_derivatives, &
                           PETSC_FALSE)
  endif

  if (hydrate_analytical_derivatives) then
    J = jac
  else
    if (hydrate_central_diff_jacobian) then
      do idof = 1, option%nflowdof
        call HydrateAccumulation(hyd_auxvar(idof),global_auxvar, &
                               material_auxvar,z,offset,meth, &
                               soil_heat_capacity,option, &
                               res_pert_plus,jac_pert,PETSC_FALSE,PETSC_FALSE)
      
        call HydrateAccumulation(hyd_auxvar(idof+option%nflowdof), &
                               global_auxvar,material_auxvar,z,offset,meth, &
                               soil_heat_capacity,option, &
                               res_pert_minus,jac_pert,PETSC_FALSE,PETSC_FALSE)
      
        do irow = 1, option%nflowdof
          J(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0 * &
                        hyd_auxvar(idof)%pert)
        enddo !irow
      enddo ! idof
    else
      do idof = 1, option%nflowdof
        call HydrateAccumulation(hyd_auxvar(idof),global_auxvar, &
                               material_auxvar,z,offset,meth, &
                               soil_heat_capacity,option, &
                               res_pert_plus,jac_pert,PETSC_FALSE,PETSC_FALSE)

        do irow = 1, option%nflowdof
          J(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                        hyd_auxvar(idof)%pert
        enddo !irow
      enddo ! idof
    endif
  endif

end subroutine HydrateAccumDerivative

! ************************************************************************** !

subroutine HydrateFluxDerivative(hyd_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 thermal_conductivity_up, &
                                 hyd_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 thermal_conductivity_dn, &
                                 area, dist, upwind_direction_, &
                                 methanogenesis, &
                                 hydrate_parameter, &
                                 option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Option_module
  use Material_Aux_class
  use Upwind_Direction_module, only : count_upwind_direction_flip
  
  implicit none
  
  type(hydrate_auxvar_type) :: hyd_auxvar_up(0:), hyd_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(methanogenesis_type) :: methanogenesis
  type(hydrate_parameter_type) :: hydrate_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_up(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_dn(option%nflowdof,option%nflowdof)
  PetscReal :: Jdummy(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
!geh:print *, 'HydrateFluxDerivative'
  option%iflag = -2

  if (.not. hydrate_central_diff_jacobian) then
    call HydrateFlux(hyd_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   thermal_conductivity_up, &
                   hyd_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   thermal_conductivity_dn, &
                   area,dist,upwind_direction_, &
                   methanogenesis, &
                   hydrate_parameter, &
                   option,v_darcy,res,Janal_up,Janal_dn,&
                   hydrate_analytical_derivatives, &
                   PETSC_FALSE, & ! update the upwind direction
                   ! avoid double counting upwind direction flip
                   PETSC_FALSE, & ! count upwind direction flip
                   PETSC_FALSE)
  endif
  if (hydrate_analytical_derivatives) then
    Jup = Janal_up
    Jdn = Janal_dn
  else
    ! upgradient derivatives
    if(hydrate_central_diff_jacobian) then
       do idof = 1, option%nflowdof
         call HydrateFlux(hyd_auxvar_up(idof),global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,res_pert_plus,Jdummy,Jdummy, &
                       PETSC_FALSE, & ! analytical derivatives
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, &
                       PETSC_FALSE)
     
         call HydrateFlux(hyd_auxvar_up(idof+option%nflowdof), &
                       global_auxvar_up,material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,res_pert_minus,Jdummy,Jdummy, &
                       PETSC_FALSE, & ! analytical derivatives
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, &
                       PETSC_FALSE) 
         do irow = 1, option%nflowdof
           Jup(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/(2.d0 * &
                            hyd_auxvar_up(idof)%pert)
         enddo !irow
       enddo ! idof
    else
      do idof = 1, option%nflowdof
        call HydrateFlux(hyd_auxvar_up(idof),global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,res_pert_plus,Jdummy,Jdummy, &
                       PETSC_FALSE, & ! analytical derivatives
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, &
                       PETSC_FALSE)

        do irow = 1, option%nflowdof
          Jup(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                            hyd_auxvar_up(idof)%pert
  !geh:print *, 'up: ', irow, idof, Jup(irow,idof), hyd_auxvar_up(idof)%pert
        enddo !irow
      enddo ! idof
    endif

    ! downgradient derivatives
    if (hydrate_central_diff_jacobian) then
      do idof = 1, option%nflowdof
        call HydrateFlux(hyd_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,res_pert_plus,Jdummy,Jdummy, &
                       PETSC_FALSE, & ! analytical derivatives
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, &
                       PETSC_FALSE)
     
        call HydrateFlux(hyd_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn(idof+option%nflowdof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,res_pert_minus,Jdummy,Jdummy, &
                       PETSC_FALSE, & ! analytical derivatives
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, &
                       PETSC_FALSE) 
      
        do irow = 1, option%nflowdof
          Jdn(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0*  &
                            hyd_auxvar_dn(idof)%pert)
  !geh:print *, 'dn: ', irow, idof, Jdn(irow,idof), hyd_auxvar_dn(idof)%pert
        enddo !irow
      enddo ! idof
    else
      do idof = 1, option%nflowdof
        call HydrateFlux(hyd_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                       material_auxvar_up, &
                       thermal_conductivity_up, &
                       hyd_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       thermal_conductivity_dn, &
                       area,dist,upwind_direction_, &
                       methanogenesis, &
                       hydrate_parameter, &
                       option,v_darcy,res_pert_plus,Jdummy,Jdummy, &
                       PETSC_FALSE, & ! analytical derivatives
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, &
                       PETSC_FALSE)

        do irow = 1, option%nflowdof
          Jdn(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                            hyd_auxvar_dn(idof)%pert
  !geh:print *, 'dn: ', irow, idof, Jdn(irow,idof), hyd_auxvar_dn(idof)%pert
        enddo !irow
      enddo ! idof
    endif
  endif
end subroutine HydrateFluxDerivative

! ************************************************************************** !

subroutine HydrateBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   hyd_auxvar_up, &
                                   global_auxvar_up, &
                                   hyd_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   thermal_conductivity_dn, &
                                   area,dist,upwind_direction_, &
                                   methanogenesis, &
                                   hydrate_parameter, &
                                   option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module 
  use Material_Aux_class
  use Upwind_Direction_module, only : count_upwind_direction_flip
  
  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(HYDRATE_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(hydrate_auxvar_type) :: hyd_auxvar_up, hyd_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(methanogenesis_type) :: methanogenesis
  type(hydrate_parameter_type) :: hydrate_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscInt :: idof, irow
  PetscReal :: Jdum(option%nflowdof,option%nflowdof)

  Jdn = 0.d0
!geh:print *, 'HydrateBCFluxDerivative'

  option%iflag = -2

  if (.not. hydrate_central_diff_jacobian) then
    call HydrateBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     hyd_auxvar_up,global_auxvar_up, &
                     hyd_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     thermal_conductivity_dn, &
                     area,dist,upwind_direction_, &
                     methanogenesis, &
                     hydrate_parameter, &
                     option,v_darcy,res,Jdum, &
                     hydrate_analytical_derivatives, &
                     PETSC_FALSE, & ! update the upwind direction
                     ! avoid double counting upwind direction flip
                     PETSC_FALSE, & ! count upwind direction flip
                     PETSC_FALSE)
  endif

  if (hydrate_analytical_derivatives) then
    Jdn = Jdum
  else
    ! downgradient derivatives
    if (hydrate_central_diff_jacobian) then
      do idof = 1, option%nflowdof
        call HydrateBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         hyd_auxvar_up,global_auxvar_up, &
                         hyd_auxvar_dn(idof),global_auxvar_dn, &
                         material_auxvar_dn, &
                         thermal_conductivity_dn, &
                         area,dist,upwind_direction_, &
                         methanogenesis, &
                         hydrate_parameter, &
                         option,v_darcy,res_pert_plus,Jdum, &
                         PETSC_FALSE, & ! analytical derivatives
                         PETSC_FALSE, & ! update the upwind direction
                         count_upwind_direction_flip, &
                         PETSC_FALSE)

        call HydrateBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         hyd_auxvar_up,global_auxvar_up, &
                         hyd_auxvar_dn(idof+option%nflowdof),global_auxvar_dn, &
                         material_auxvar_dn, &
                         thermal_conductivity_dn, &
                         area,dist,upwind_direction_, &
                         methanogenesis, &
                         hydrate_parameter, &
                         option,v_darcy,res_pert_minus,Jdum, &
                         PETSC_FALSE, & ! analytical derivatives
                         PETSC_FALSE, & ! update the upwind direction
                         count_upwind_direction_flip, &
                         PETSC_FALSE)
    
      
        do irow = 1, option%nflowdof
          Jdn(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0 * &
                            hyd_auxvar_dn(idof)%pert)
        enddo !irow
      enddo ! idof
    else
      do idof = 1, option%nflowdof
        call HydrateBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         hyd_auxvar_up,global_auxvar_up, &
                         hyd_auxvar_dn(idof),global_auxvar_dn, &
                         material_auxvar_dn, &
                         thermal_conductivity_dn, &
                         area,dist,upwind_direction_, &
                         methanogenesis, &
                         hydrate_parameter, &
                         option,v_darcy,res_pert_plus,Jdum, &
                         PETSC_FALSE, & ! analytical derivatives
                         PETSC_FALSE, & ! update the upwind direction
                         count_upwind_direction_flip, &
                         PETSC_FALSE)

        do irow = 1, option%nflowdof
          Jdn(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                            hyd_auxvar_dn(idof)%pert
        enddo !irow
      enddo ! idof
    endif
  endif

end subroutine HydrateBCFluxDerivative

! ************************************************************************** !

subroutine HydrateSrcSinkDerivative(option,source_sink,hyd_auxvar_ss, &
                                    hyd_auxvars,global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module
  use Coupler_module

  implicit none

  type(option_type) :: option
  type(coupler_type), pointer :: source_sink
  type(hydrate_auxvar_type) :: hyd_auxvars(0:), hyd_auxvar_ss
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: qsrc(3)
  PetscInt :: flow_src_sink_type
  PetscReal :: res(option%nflowdof), res_pert_plus(option%nflowdof)
  PetscReal :: res_pert_minus(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow
  PetscReal :: Jdum(option%nflowdof,option%nflowdof)
  
  qsrc = source_sink%flow_condition%hydrate%rate%dataset%rarray(:)
  flow_src_sink_type = source_sink%flow_condition%hydrate%rate%itype

  option%iflag = -3
  
  if (.not. hydrate_central_diff_jacobian) then
    call HydrateSrcSink(option,qsrc,flow_src_sink_type,hyd_auxvar_ss, &
                      hyd_auxvars(ZERO_INTEGER),global_auxvar,dummy_real, &
                      scale,res,Jdum,hydrate_analytical_derivatives, &
                      PETSC_FALSE)
  endif

  if (hydrate_analytical_derivatives) then
    Jac = Jdum
  else                      
    ! downgradient derivatives
    if (hydrate_central_diff_jacobian) then
      do idof = 1, option%nflowdof
        call HydrateSrcSink(option,qsrc,flow_src_sink_type,hyd_auxvar_ss, &
                          hyd_auxvars(idof),global_auxvar,dummy_real, &
                          scale,res_pert_plus,Jdum,PETSC_FALSE,PETSC_FALSE)
        call HydrateSrcSink(option,qsrc,flow_src_sink_type,hyd_auxvar_ss, &
                          hyd_auxvars(idof+option%nflowdof),global_auxvar, &
                          dummy_real,scale,res_pert_minus,Jdum,PETSC_FALSE, &
                          PETSC_FALSE)
      
        do irow = 1, option%nflowdof
          Jac(irow,idof) = (res_pert_plus(irow)-res_pert_minus(irow))/ (2.d0 * &
                          hyd_auxvars(idof)%pert)
        enddo !irow
      enddo ! idof
    else
      do idof = 1, option%nflowdof
        call HydrateSrcSink(option,qsrc,flow_src_sink_type,hyd_auxvar_ss, &
                          hyd_auxvars(idof),global_auxvar,dummy_real, &
                          scale,res_pert_plus,Jdum,PETSC_FALSE,PETSC_FALSE)

        do irow = 1, option%nflowdof
          Jac(irow,idof) = (res_pert_plus(irow)-res(irow))/ &
                          hyd_auxvars(idof)%pert
        enddo !irow
      enddo ! idof 
    endif
  endif
  
end subroutine HydrateSrcSinkDerivative

! ************************************************************************** !

function HydrateAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: HydrateAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  if (iphase == LIQUID_PHASE) then
    if (istate_up == G_STATE) then
      HydrateAverageDensity = density_dn(iphase)
      dden_dn = 1.d0
    else if (istate_dn == G_STATE) then
      HydrateAverageDensity = density_up(iphase)
      dden_up = 1.d0
    else
      HydrateAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == L_STATE) then
      HydrateAverageDensity = density_dn(iphase)
      dden_dn = 1.d0      
    else if (istate_dn == L_STATE) then
      HydrateAverageDensity = density_up(iphase)
      dden_up = 1.d0      
    else
      HydrateAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0      
    endif
  endif

end function HydrateAverageDensity

! ************************************************************************** !

subroutine HydrateAuxVarDiff(idof,hydrate_auxvar,global_auxvar, &
                             material_auxvar, &
                             hydrate_auxvar_pert,global_auxvar_pert, &
                             material_auxvar_pert, &
                             pert,string,compare_analytical_derivative, &
                             option)

  use Option_module
  use Global_Aux_module
  use Material_Aux_class  

  implicit none
  
  type(option_type) :: option
  PetscInt :: idof
  type(hydrate_auxvar_type) :: hydrate_auxvar, hydrate_auxvar_pert
  type(global_auxvar_type) :: global_auxvar, global_auxvar_pert
  class(material_auxvar_type) :: material_auxvar, material_auxvar_pert
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: compare_analytical_derivative
  PetscReal :: pert

  
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation
  PetscReal :: liquid_mass_pert, gas_mass_pert
  PetscReal :: liquid_density_pert, gas_density_pert
  PetscReal :: liquid_energy_pert, gas_energy_pert
  PetscReal :: liquid_saturation_pert, gas_saturation_pert
  
  PetscReal :: dpl 
  PetscReal :: dpg 
  PetscReal :: dpa 
  PetscReal :: dpc 
  PetscReal :: dpv 
  PetscReal :: dps 
  PetscReal :: dsatl
  PetscReal :: dsatg
  PetscReal :: ddenl  
  PetscReal :: ddeng  
  PetscReal :: ddenlkg
  PetscReal :: ddengkg
  PetscReal :: dUl 
  PetscReal :: dHl  
  PetscReal :: dUg  
  PetscReal :: dHg  
  PetscReal :: dUv
  PetscReal :: dHv  
  PetscReal :: dUa  
  PetscReal :: dHa  
  PetscReal :: dpsat  
  PetscReal :: dmobilityl  
  PetscReal :: dmobilityg  
  PetscReal :: dxmolwl
  PetscReal :: dxmolal
  PetscReal :: dxmolwg
  PetscReal :: dxmolag
  PetscReal :: denv
  PetscReal :: dena
  PetscReal :: dHc
  PetscReal :: dmug
  
  PetscReal, parameter :: uninitialized_value = -999.d0
  
  dpl = uninitialized_value
  dpg = uninitialized_value
  dpa = uninitialized_value
  dpc = uninitialized_value
  dpv = uninitialized_value
  dps = uninitialized_value
  dsatl = uninitialized_value
  dsatg = uninitialized_value
  ddenl = uninitialized_value
  ddeng = uninitialized_value
  ddenlkg = uninitialized_value
  ddengkg = uninitialized_value
  dUl = uninitialized_value
  dHl = uninitialized_value
  dUg = uninitialized_value
  dHg = uninitialized_value
  dUv = uninitialized_value
  dHv = uninitialized_value
  dUa = uninitialized_value
  dHa = uninitialized_value
  dpsat = uninitialized_value
  dmobilityl = uninitialized_value
  dmobilityg = uninitialized_value
  dxmolwl = uninitialized_value
  dxmolal = uninitialized_value
  dxmolwg = uninitialized_value
  dxmolag = uninitialized_value
  denv = uninitialized_value
  dena = uninitialized_value
  dHc = uninitialized_value
  dmug = uninitialized_value

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0
    
  if (compare_analytical_derivative) then
    select case(global_auxvar%istate)
      case(L_STATE)
        select case(idof)
          case(1) ! liquid pressure pl
            dpl = 1.d0
            dpv = 0.d0
            dpa = hydrate_auxvar%d%Hc_p*hydrate_auxvar%xmol(acid,lid)
            dps = 0.d0
            dpv = hydrate_auxvar%d%pv_p
            ddenl = hydrate_auxvar%d%denl_pl
            dUl = hydrate_auxvar%d%Ul_pl
            dHl = hydrate_auxvar%d%Hl_pl
            dmobilityl = hydrate_auxvar%d%mobilityl_pl
          case(2) ! xmole air in liquid
            dxmolwl = -1.d0
            dxmolal = 1.d0
          case(3) ! temperature
            dpl = 0.d0 ! pl = pg - pc
            dpg = 0.d0 ! pg = pg
            dpa = hydrate_auxvar%d%Hc_T*hydrate_auxvar%xmol(acid,lid)
            dps = hydrate_auxvar%d%psat_T
            dHc = hydrate_auxvar%d%Hc_T            
            dsatl = 0.d0
            dsatg = 0.d0            
            ddenl = hydrate_auxvar%d%denl_T
            ddeng = hydrate_auxvar%d%deng_T
            ddenlkg = ddenl*fmw_comp(1)
            ddengkg = hydrate_auxvar%d%dengkg_T
            dUl = hydrate_auxvar%d%Ul_T
            dHl = hydrate_auxvar%d%Hl_T
            
            dpv = hydrate_auxvar%d%pv_T
            dmobilityl = hydrate_auxvar%d%mobilityl_T
            dxmolwl = hydrate_auxvar%d%xmol_T(wid,lid)
            dxmolal = hydrate_auxvar%d%xmol_T(acid,lid)
        end select
      case(G_STATE)
        select case(idof)
          case(1)
            dpg = 1.d0 ! pg = pg
            dpv = hydrate_auxvar%d%pv_p
            dps = 0.d0
            dHc = hydrate_auxvar%d%Hc_p
            ddeng = hydrate_auxvar%d%deng_pg
            ddengkg = hydrate_auxvar%d%dengkg_pg
            dUg = hydrate_auxvar%d%Ug_pg
            dHg = hydrate_auxvar%d%Hg_pg
            
            dHv = hydrate_auxvar%d%Hv_pg
            dUv = hydrate_auxvar%d%Uv_pg
            dHa = hydrate_auxvar%d%Ha_pg
            dUa = hydrate_auxvar%d%Ua_pg
            
            dmug = hydrate_auxvar%d%mug_pg
            dmobilityg = hydrate_auxvar%d%mobilityg_pg
            dxmolwg = hydrate_auxvar%d%xmol_p(wid,gid)
            dxmolag = hydrate_auxvar%d%xmol_p(acid,gid)          
          case(2)
            dpg = 0.d0
            dpa = 1.d0
            ddeng = hydrate_auxvar%d%deng_pa
            dpv = hydrate_auxvar%d%pv_pa
            dUg = hydrate_auxvar%d%Ug_pa
            dHg = hydrate_auxvar%d%Hg_pa
            dHv = hydrate_auxvar%d%Hv_pa
            dUv = hydrate_auxvar%d%Uv_pa
            dHa = hydrate_auxvar%d%Ha_pa
            dUa = hydrate_auxvar%d%Ua_pa
            ! for gas state, derivative wrt air pressure is under lid
            dxmolwg = hydrate_auxvar%d%xmol_p(wid,lid)
            dxmolag = hydrate_auxvar%d%xmol_p(acid,lid)          
            dmobilityg = hydrate_auxvar%d%mobilityg_pa
          case(3)
            dpg = 0.d0
            dpa = 0.d0
            dpv = 0.d0
            dps = hydrate_auxvar%d%psat_T
            dHc = hydrate_auxvar%d%Hc_T            
            ddeng = hydrate_auxvar%d%deng_T
            ddengkg = hydrate_auxvar%d%dengkg_T
            dUg = hydrate_auxvar%d%Ug_T
            dHg = hydrate_auxvar%d%Hg_T
            
            dHv = hydrate_auxvar%d%Hv_T
            dUv = hydrate_auxvar%d%Uv_T
            dHa = hydrate_auxvar%d%Ha_T
            dUa = hydrate_auxvar%d%Ua_T
            denv = hydrate_auxvar%d%denv_T
            dena = hydrate_auxvar%d%dena_T
            
            dmug = hydrate_auxvar%d%mug_T
            dmobilityg = hydrate_auxvar%d%mobilityg_T
            dxmolwg = hydrate_auxvar%d%xmol_T(wid,gid)
            dxmolag = hydrate_auxvar%d%xmol_T(acid,gid)          
        end select
      case(GA_STATE)
        select case(idof)
          case(1) ! gas pressure pg
            dpl = 1.d0 ! pl = pg - pc
            dpg = 1.d0 ! pg = pg
            dpa = 1.d0 ! pa = pg - pv
            dpv = 0.d0
            dps = 0.d0
            dsatl = 0.d0
            dsatg = 0.d0
            dHc = hydrate_auxvar%d%Hc_p
            ddenl = hydrate_auxvar%d%denl_pl*dpl
            ddeng = hydrate_auxvar%d%deng_pg
            ddenlkg = ddenl*fmw_comp(1)
            ddengkg = hydrate_auxvar%d%dengkg_pg
            dUl = hydrate_auxvar%d%Ul_pl
            dHl = hydrate_auxvar%d%Hl_pl
            dUg = hydrate_auxvar%d%Ug_pg
            dHg = hydrate_auxvar%d%Hg_pg
            
            denv = hydrate_auxvar%d%denv_pg
            dena = hydrate_auxvar%d%dena_pg
            dHv = hydrate_auxvar%d%Hv_pg
            dUv = hydrate_auxvar%d%Uv_pg
            dHa = hydrate_auxvar%d%Ha_pg
            dUa = hydrate_auxvar%d%Ua_pg
            
            dmug = hydrate_auxvar%d%mug_pg
            dmobilityl = hydrate_auxvar%d%mobilityl_pl
            dmobilityg = hydrate_auxvar%d%mobilityg_pg
            dxmolwl = hydrate_auxvar%d%xmol_p(wid,lid)
            dxmolal = hydrate_auxvar%d%xmol_p(acid,lid)
            dxmolwg = hydrate_auxvar%d%xmol_p(wid,gid)
            dxmolag = hydrate_auxvar%d%xmol_p(acid,gid)
          case(2) ! gas saturation
            dpl = -1.d0*hydrate_auxvar%d%pc_satg ! pl = pg - pc
            dpg = 0.d0
            dpa = 0.d0
            dpc = hydrate_auxvar%d%pc_satg
            dpv = 0.d0 
            dps = 0.d0
            dsatl = -1.d0
            dsatg = 1.d0
            ddenl = 0.d0
            dmobilityl = hydrate_auxvar%d%mobilityl_satg
            dmobilityg = hydrate_auxvar%d%mobilityg_satg
          case(3) ! temperature
            dpl = 0.d0 ! pl = pg - pc
            dpg = 0.d0 ! pg = pg
            dpa = -1.d0*hydrate_auxvar%d%psat_T ! pa = pg - pv
            dpv = hydrate_auxvar%d%psat_T
            dps = hydrate_auxvar%d%psat_T
            dHc = hydrate_auxvar%d%Hc_T            
            dsatl = 0.d0
            dsatg = 0.d0            
            ddenl = hydrate_auxvar%d%denl_T
            ddeng = hydrate_auxvar%d%deng_T
            ddenlkg = ddenl*fmw_comp(1)
            ddengkg = hydrate_auxvar%d%dengkg_T
            dUl = hydrate_auxvar%d%Ul_T
            dHl = hydrate_auxvar%d%Hl_T
            dUg = hydrate_auxvar%d%Ug_T
            dHg = hydrate_auxvar%d%Hg_T
            
            dHv = hydrate_auxvar%d%Hv_T
            dUv = hydrate_auxvar%d%Uv_T
            dHa = hydrate_auxvar%d%Ha_T
            dUa = hydrate_auxvar%d%Ua_T
            denv = hydrate_auxvar%d%denv_T
            dena = hydrate_auxvar%d%dena_T
            
            dmug = hydrate_auxvar%d%mug_T
            dmobilityl = hydrate_auxvar%d%mobilityl_T
            dmobilityg = hydrate_auxvar%d%mobilityg_T
            dxmolwl = hydrate_auxvar%d%xmol_T(wid,lid)
            dxmolal = hydrate_auxvar%d%xmol_T(acid,lid)
            dxmolwg = hydrate_auxvar%d%xmol_T(wid,gid)
            dxmolag = hydrate_auxvar%d%xmol_T(acid,gid)
        end select
      end select
    endif

  print *, '--------------------------------------------------------'
  print *, 'Derivative with respect to ' // trim(string)
  select case(global_auxvar%istate)
    case(L_STATE)
      print *, '     Thermodynamic state: Liquid phase'
      liquid_density = hydrate_auxvar%den(lid)
      liquid_energy = hydrate_auxvar%U(lid)
      liquid_saturation = hydrate_auxvar%sat(lid)
    case(G_STATE)
      print *, '     Thermodynamic state: Gas phase'
      gas_density = hydrate_auxvar%den(gid)
      gas_energy = hydrate_auxvar%U(gid)
      gas_saturation = hydrate_auxvar%sat(gid)
    case(GA_STATE)
      print *, '     Thermodynamic state: Two phase'
      liquid_density = hydrate_auxvar%den(lid)
      gas_density = hydrate_auxvar%den(gid)
      liquid_energy = hydrate_auxvar%U(lid)
      gas_energy = hydrate_auxvar%U(gid)
      liquid_saturation = hydrate_auxvar%sat(lid)
      gas_saturation = hydrate_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*hydrate_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*hydrate_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 hydrate_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*hydrate_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*hydrate_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              hydrate_auxvar%effective_porosity*material_auxvar%volume
  select case(global_auxvar_pert%istate)
    case(L_STATE)
      print *, '     Thermodynamic state (pert): Liquid phase'
      liquid_density_pert = hydrate_auxvar_pert%den(lid)
      liquid_energy_pert = hydrate_auxvar_pert%U(lid)
      liquid_saturation_pert = hydrate_auxvar_pert%sat(lid)
      gas_density_pert = 0.d0
      gas_energy_pert = 0.d0
      gas_saturation_pert = 0.d0
    case(G_STATE)
      print *, '     Thermodynamic state (pert): Gas phase'
      liquid_density_pert = 0.d0
      liquid_energy_pert = 0.d0
      liquid_saturation_pert = 0.d0
      gas_density_pert = hydrate_auxvar_pert%den(gid)
      gas_energy_pert = hydrate_auxvar_pert%U(gid)
      gas_saturation_pert = hydrate_auxvar_pert%sat(gid)
    case(GA_STATE)
      print *, '     Thermodynamic state (pert): Two phase'
      liquid_density_pert = hydrate_auxvar_pert%den(lid)
      gas_density_pert = hydrate_auxvar_pert%den(gid)
      liquid_energy_pert = hydrate_auxvar_pert%U(lid)
      gas_energy_pert = hydrate_auxvar_pert%U(gid)
      liquid_saturation_pert = hydrate_auxvar_pert%sat(lid)
      gas_saturation_pert = hydrate_auxvar_pert%sat(gid)
  end select  
  liquid_mass_pert = (liquid_density_pert*hydrate_auxvar_pert%xmol(lid,lid)* & 
                 liquid_saturation_pert+ &
                 gas_density_pert*hydrate_auxvar_pert%xmol(lid,gid)* & 
                 gas_saturation_pert)* & 
                 hydrate_auxvar_pert%effective_porosity*material_auxvar_pert%volume
  gas_mass_pert = (liquid_density_pert*hydrate_auxvar_pert%xmol(gid,lid)* & 
              liquid_saturation_pert+ &
              gas_density_pert*hydrate_auxvar_pert%xmol(gid,gid)* & 
              gas_saturation_pert)* & 
              hydrate_auxvar_pert%effective_porosity*material_auxvar_pert%volume 
              
              
  call HydrateAuxVarPrintResult('tot liq comp mass [kmol]', &
                                (liquid_mass_pert-liquid_mass)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call HydrateAuxVarPrintResult('tot gas comp mass [kmol]', &
                                (gas_mass_pert-gas_mass)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call HydrateAuxVarPrintResult('             energy [MJ]', &
                                ((liquid_mass_pert*liquid_energy_pert + &
                                  gas_mass_pert*gas_energy_pert)- &
                                 (liquid_mass*liquid_energy + &
                                  gas_mass*gas_energy))/pert, &
                                uninitialized_value,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         liquid pressure', &
                                (hydrate_auxvar_pert%pres(lid)-hydrate_auxvar%pres(lid))/pert, &
                                dpl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('            gas pressure', &
                                (hydrate_auxvar_pert%pres(gid)-hydrate_auxvar%pres(gid))/pert, &
                                dpg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('            air pressure', &
                                (hydrate_auxvar_pert%pres(apid)-hydrate_auxvar%pres(apid))/pert, &
                                dpa,uninitialized_value,option)
  call HydrateAuxVarPrintResult('      capillary pressure', &
                                (hydrate_auxvar_pert%pres(cpid)-hydrate_auxvar%pres(cpid))/pert, &
                                dpc,uninitialized_value,option)
  call HydrateAuxVarPrintResult('          vapor pressure', &
                                (hydrate_auxvar_pert%pres(vpid)-hydrate_auxvar%pres(vpid))/pert, &
                                dpv,uninitialized_value,option)
  call HydrateAuxVarPrintResult("        Henry's constant", &
                                (hydrate_auxvar_pert%d%Hc-hydrate_auxvar%d%Hc)/pert, &
                                dHc,uninitialized_value,option)
  call HydrateAuxVarPrintResult('     saturation pressure', &
                                (hydrate_auxvar_pert%pres(spid)-hydrate_auxvar%pres(spid))/pert, &
                                dps,uninitialized_value,option)
  call HydrateAuxVarPrintResult('       liquid saturation', &
                                (hydrate_auxvar_pert%sat(lid)-hydrate_auxvar%sat(lid))/pert, &
                                dsatl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('          gas saturation', &
                                (hydrate_auxvar_pert%sat(gid)-hydrate_auxvar%sat(gid))/pert, &
                                dsatg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('   liquid density [kmol]', &
                                (hydrate_auxvar_pert%den(lid)-hydrate_auxvar%den(lid))/pert, &
                                ddenl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('      gas density [kmol]', &
                                (hydrate_auxvar_pert%den(gid)-hydrate_auxvar%den(gid))/pert, &
                                ddeng,uninitialized_value,option)
  call HydrateAuxVarPrintResult('     liquid density [kg]', &
                                (hydrate_auxvar_pert%den_kg(lid)-hydrate_auxvar%den_kg(lid))/pert, &
                                ddenlkg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('        gas density [kg]', &
                                (hydrate_auxvar_pert%den_kg(gid)-hydrate_auxvar%den_kg(gid))/pert, &
                                ddengkg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         temperature [C]', &
                                (hydrate_auxvar_pert%temp-hydrate_auxvar%temp)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call HydrateAuxVarPrintResult('      liquid H [MJ/kmol]', &
                                (hydrate_auxvar_pert%H(lid)-hydrate_auxvar%H(lid))/pert, &
                                dHl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         gas H [MJ/kmol]', &
                                (hydrate_auxvar_pert%H(gid)-hydrate_auxvar%H(gid))/pert, &
                                dHg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('      liquid U [MJ/kmol]', &
                                (hydrate_auxvar_pert%U(lid)-hydrate_auxvar%U(lid))/pert, &
                                dUl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         gas U [MJ/kmol]', &
                                (hydrate_auxvar_pert%U(gid)-hydrate_auxvar%U(gid))/pert, &
                                dUg,uninitialized_value,option)
  !------------------------------
  call HydrateAuxVarPrintResult('       vapor H [MJ/kmol]', &
                                (hydrate_auxvar_pert%d%Hv-hydrate_auxvar%d%Hv)/pert, &
                                dHv,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         air H [MJ/kmol]', &
                                (hydrate_auxvar_pert%d%Ha-hydrate_auxvar%d%Ha)/pert, &
                                dHa,uninitialized_value,option)
  call HydrateAuxVarPrintResult('       vapor U [MJ/kmol]', &
                                (hydrate_auxvar_pert%d%Uv-hydrate_auxvar%d%Uv)/pert, &
                                dUv,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         air U [MJ/kmol]', &
                                (hydrate_auxvar_pert%d%Ua-hydrate_auxvar%d%Ua)/pert, &
                                dUa,uninitialized_value,option)
  call HydrateAuxVarPrintResult('    vapor density [kmol]', &
                                (hydrate_auxvar_pert%d%denv-hydrate_auxvar%d%denv)/pert, &
                                denv,uninitialized_value,option)
  call HydrateAuxVarPrintResult('      air density [kmol]', &
                                (hydrate_auxvar_pert%d%dena-hydrate_auxvar%d%dena)/pert, &
                                dena,uninitialized_value,option)
  !------------------------------                                
  call HydrateAuxVarPrintResult('     X (water in liquid)', &
                                (hydrate_auxvar_pert%xmol(wid,lid)-hydrate_auxvar%xmol(wid,lid))/pert, &
                                dxmolwl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('       X (air in liquid)', &
                                (hydrate_auxvar_pert%xmol(acid,lid)-hydrate_auxvar%xmol(acid,lid))/pert, &
                                dxmolal,uninitialized_value,option)
  call HydrateAuxVarPrintResult('        X (water in gas)', &
                                (hydrate_auxvar_pert%xmol(wid,gid)-hydrate_auxvar%xmol(wid,gid))/pert, &
                                dxmolwg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('          X (air in gas)', &
                                (hydrate_auxvar_pert%xmol(acid,gid)-hydrate_auxvar%xmol(acid,gid))/pert, &
                                dxmolag,uninitialized_value,option)
  call HydrateAuxVarPrintResult('         liquid mobility', &
                                (hydrate_auxvar_pert%mobility(lid)-hydrate_auxvar%mobility(lid))/pert, &
                                dmobilityl,uninitialized_value,option)
  call HydrateAuxVarPrintResult('            gas mobility', &
                                (hydrate_auxvar_pert%mobility(gid)-hydrate_auxvar%mobility(gid))/pert, &
                                dmobilityg,uninitialized_value,option)
  call HydrateAuxVarPrintResult('           gas viscosity', &
                                (hydrate_auxvar_pert%d%mug-hydrate_auxvar%d%mug)/pert, &
                                dmug,uninitialized_value,option)
  call HydrateAuxVarPrintResult('      effective porosity', &
                                (hydrate_auxvar_pert%effective_porosity-hydrate_auxvar%effective_porosity)/pert, &
                                uninitialized_value,uninitialized_value,option)
#if 0                                
100 format(a,2(es13.5),es16.8)  
  write(*,100) 'tot liq comp mass [kmol]: ', (liquid_mass_pert-liquid_mass)/pert
  write(*,100) 'tot gas comp mass [kmol]: ', (gas_mass_pert-gas_mass)/pert
  write(*,100) '             energy [MJ]: ', ((liquid_mass_pert*liquid_energy_pert + &
                                           gas_mass_pert*gas_energy_pert)- &
                                          (liquid_mass*liquid_energy + &
                                           gas_mass*gas_energy))/pert
  write(*,100) '         liquid pressure: ', (hydrate_auxvar_pert%pres(lid)-hydrate_auxvar%pres(lid))/pert,dpl
  write(*,100) '            gas pressure: ', (hydrate_auxvar_pert%pres(gid)-hydrate_auxvar%pres(gid))/pert,dpg
  write(*,100) '            air pressure: ', (hydrate_auxvar_pert%pres(apid)-hydrate_auxvar%pres(apid))/pert,dpa !,hydrate_auxvar_pert%pres(apid)-hydrate_auxvar%pres(apid)
  write(*,100) '      capillary pressure: ', (hydrate_auxvar_pert%pres(cpid)-hydrate_auxvar%pres(cpid))/pert,dpc
  write(*,100) '          vapor pressure: ', (hydrate_auxvar_pert%pres(vpid)-hydrate_auxvar%pres(vpid))/pert,dpv !,hydrate_auxvar_pert%pres(vpid)-hydrate_auxvar%pres(vpid)
  write(*,100) "        Henry's constant: ", (hydrate_auxvar_pert%d%Hc-hydrate_auxvar%d%Hc)/pert,dHc
  write(*,100) '     saturation pressure: ', (hydrate_auxvar_pert%pres(spid)-hydrate_auxvar%pres(spid))/pert,dps
  write(*,100) '       liquid saturation: ', (hydrate_auxvar_pert%sat(lid)-hydrate_auxvar%sat(lid))/pert,dsatl
  write(*,100) '          gas saturation: ', (hydrate_auxvar_pert%sat(gid)-hydrate_auxvar%sat(gid))/pert,dsatg
  write(*,100) '   liquid density [kmol]: ', (hydrate_auxvar_pert%den(lid)-hydrate_auxvar%den(lid))/pert,ddenl
  write(*,100) '      gas density [kmol]: ', (hydrate_auxvar_pert%den(gid)-hydrate_auxvar%den(gid))/pert,ddeng
  write(*,100) '     liquid density [kg]: ', (hydrate_auxvar_pert%den_kg(lid)-hydrate_auxvar%den_kg(lid))/pert,ddenl*fmw_comp(1)
  write(*,100) '        gas density [kg]: ', (hydrate_auxvar_pert%den_kg(gid)-hydrate_auxvar%den_kg(gid))/pert,ddengkg
  write(*,100) '         temperature [C]: ', (hydrate_auxvar_pert%temp-hydrate_auxvar%temp)/pert
  write(*,100) '      liquid H [MJ/kmol]: ', (hydrate_auxvar_pert%H(lid)-hydrate_auxvar%H(lid))/pert,dHl
  write(*,100) '         gas H [MJ/kmol]: ', (hydrate_auxvar_pert%H(gid)-hydrate_auxvar%H(gid))/pert,dHg
  write(*,100) '      liquid U [MJ/kmol]: ', (hydrate_auxvar_pert%U(lid)-hydrate_auxvar%U(lid))/pert,dUl
  write(*,100) '         gas U [MJ/kmol]: ', (hydrate_auxvar_pert%U(gid)-hydrate_auxvar%U(gid))/pert,dUg

  write(*,100) '       vapor H [MJ/kmol]: ', (hydrate_auxvar_pert%d%Hv-hydrate_auxvar%d%Hv)/pert,dHv
  write(*,100) '         air H [MJ/kmol]: ', (hydrate_auxvar_pert%d%Ha-hydrate_auxvar%d%Ha)/pert,dHa
  write(*,100) '       vapor U [MJ/kmol]: ', (hydrate_auxvar_pert%d%Uv-hydrate_auxvar%d%Uv)/pert,dUv

  write(*,100) '         air U [MJ/kmol]: ', (hydrate_auxvar_pert%d%Ua-hydrate_auxvar%d%Ua)/pert,dUa
  write(*,100) '    vapor density [kmol]: ', (hydrate_auxvar_pert%d%denv-hydrate_auxvar%d%denv)/pert,denv
  write(*,100) '      air density [kmol]: ', (hydrate_auxvar_pert%d%dena-hydrate_auxvar%d%dena)/pert,dena

  write(*,100) '     X (water in liquid): ', (hydrate_auxvar_pert%xmol(wid,lid)-hydrate_auxvar%xmol(wid,lid))/pert,dxmolwl
  write(*,100) '       X (air in liquid): ', (hydrate_auxvar_pert%xmol(acid,lid)-hydrate_auxvar%xmol(acid,lid))/pert,dxmolal
  write(*,100) '        X (water in gas): ', (hydrate_auxvar_pert%xmol(wid,gid)-hydrate_auxvar%xmol(wid,gid))/pert,dxmolwg
  write(*,100) '          X (air in gas): ', (hydrate_auxvar_pert%xmol(acid,gid)-hydrate_auxvar%xmol(acid,gid))/pert,dxmolag
  write(*,100) '         liquid mobility: ', (hydrate_auxvar_pert%mobility(lid)-hydrate_auxvar%mobility(lid))/pert,dmobilityl
  write(*,100) '            gas mobility: ', (hydrate_auxvar_pert%mobility(gid)-hydrate_auxvar%mobility(gid))/pert,dmobilityg
  write(*,100) '      effective porosity: ', (hydrate_auxvar_pert%effective_porosity-hydrate_auxvar%effective_porosity)/pert
#endif
  write(*,*) '--------------------------------------------------------'  
  
end subroutine HydrateAuxVarDiff

! ************************************************************************** !

subroutine HydrateAuxVarPrintResult(string,numerical,analytical, &
                                    uninitialized_value,option)

  use Option_module
  use Utility_module
  
  implicit none
  
  character(len=*) :: string
  PetscReal :: numerical
  PetscReal :: analytical
  PetscReal :: uninitialized_value
  type(option_type) :: option
  
  character(len=8) :: word
  character(len=2) :: precision
  PetscReal :: tempreal
  PetscReal, parameter :: tol = 1.d-5
          
100 format(a24,': ',2(es13.5),2x,a2,x,a8,x,es16.8)

  precision = DigitsOfAccuracy(numerical,analytical)
  word = ''
  if (dabs(analytical-uninitialized_value) > 1.d-20) then
    if (dabs(analytical) > 0.d0) then
      tempreal = dabs((numerical-analytical)/analytical)
      if (tempreal < tol) then
        word = ' PASS'
      else
        word = '-FAIL-'
        if (tempreal < 1.d1*tol) word = trim(word) // ' *'
      endif
    else
      if (dabs(numerical) > 1.d-20) then
        word = '-FAIL-'
      else
        word = ' PASS'
      endif
    endif
    write(*,100) trim(string), numerical, analytical, precision, word
  else
    write(*,100) trim(string), numerical
  endif
  

end subroutine HydrateAuxVarPrintResult

! ************************************************************************** !

subroutine HydrateDiffJacobian(string,numerical_jacobian,analytical_jacobian, &
                               residual,residual_pert,perturbation, &
                               perturbation_tolerance,hydrate_auxvar,option)

  use Option_module
  use Utility_module
  
  implicit none
  
  character(len=*) :: string
  PetscReal :: numerical_jacobian(3,3)
  PetscReal :: analytical_jacobian(3,3)
  PetscReal :: residual(3)
  PetscReal :: residual_pert(3,3)
  PetscReal :: perturbation(3)
  PetscReal :: perturbation_tolerance
  type(hydrate_auxvar_type) :: hydrate_auxvar(0:)
  type(option_type) :: option
  
  PetscInt :: irow, icol
  
100 format(2i2,2es13.5,x,a2,es16.8)

  if (len_trim(string) > 1) then
    write(*,'(x,a)') string
  endif
  write(*,'(" Perturbation tolerance: ",es12.4)') perturbation_tolerance
  write(*,'(" r c  numerical    analytical   digits of accuracy")')
  do icol = 1, 3
    do irow = 1, 3
      write(*,100) irow, icol, numerical_jacobian(irow,icol), &
                   analytical_jacobian(irow,icol), &
                   DigitsOfAccuracy(numerical_jacobian(irow,icol), &
                                    analytical_jacobian(irow,icol))
    enddo
  enddo

#if 0
200 format(2es20.12)
300 format(a24,10es20.12)
  do icol = 1, 3
    write(*,'(/," dof = ",i1,"  perturbation = ",es13.5)') icol, perturbation(icol)
!    write(*,300) 'density', hydrate_auxvar(icol)%den(:), hydrate_auxvar(0)%den(:)
!    write(*,300) 'energy', hydrate_auxvar(icol)%U(:), hydrate_auxvar(0)%U(:)
    write(*,'("  residual_pert       residual")')
    do irow = 1, 3
      write(*,200) residual_pert(irow,icol), residual(irow)
    enddo
  enddo
#endif  
  
end subroutine HydrateDiffJacobian

! ************************************************************************** !

end module Hydrate_Common_module
