module WIPP_Flow_Common_module

  use WIPP_Flow_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use petscsys

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  procedure(FluxDummy), pointer :: XXFlux => null()
  procedure(BCFluxDummy), pointer :: XXBCFlux => null()

  interface
    subroutine FluxDummy(wippflo_auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         wippflo_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area, dist, upwind_direction_, &
                         wippflo_parameter, &
                         option,v_darcy,Res, &
                         derivative_call, &
                         update_upwind_direction_, &
                         count_upwind_direction_flip_, &
                         debug_connection)
      use WIPP_Flow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      implicit none
      type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy(option%nphase)
      PetscReal :: area
      PetscReal :: dist(-1:3)
      PetscInt :: upwind_direction_(option%nphase)
      type(wippflo_parameter_type) :: wippflo_parameter
      PetscReal :: Res(option%nflowdof)
      PetscBool :: derivative_call
      PetscBool :: update_upwind_direction_
      PetscBool :: count_upwind_direction_flip_
      PetscBool :: debug_connection
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                         wippflo_auxvar_up,global_auxvar_up, &
                         wippflo_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist,upwind_direction_, &
                         wippflo_parameter, &
                         option,v_darcy,Res, &
                         derivative_call, &
                         update_upwind_direction_, &
                         count_upwind_direction_flip_, &
                         debug_connection)
      use WIPP_Flow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      implicit none
      type(option_type) :: option
      PetscInt :: ibndtype(1:option%nflowdof)
      PetscInt :: auxvar_mapping(WIPPFLO_MAX_INDEX)
      PetscReal :: auxvars(:) ! from aux_real_var array
      type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_dn
      PetscReal :: area
      PetscReal :: dist(-1:3)
      PetscInt :: upwind_direction_(option%nphase)
      type(wippflo_parameter_type) :: wippflo_parameter
      PetscReal :: v_darcy(option%nphase)
      PetscReal :: Res(1:option%nflowdof)
      PetscBool :: derivative_call
      PetscBool :: update_upwind_direction_
      PetscBool :: count_upwind_direction_flip_
      PetscBool :: debug_connection
    end subroutine BCFluxDummy
  end interface

  public :: WIPPFloAccumulation, &
            WIPPFloFluxHarmonicPermOnly, &
            WIPPFloFluxLumpedHarmonic, &
            WIPPFloBCFluxHarmonicPermOnly, &
            WIPPFloBCFluxLumpedHarmonic, &
            WIPPFloSrcSink, &
            WIPPFloAccumDerivative, &
            XXFluxDerivative, &
            XXBCFluxDerivative, &
            WIPPFloSrcSinkDerivative, &
            WIPPFloAverageDensity

  public :: XXFlux, &
            XXBCFlux
            
contains

! ************************************************************************** !

subroutine WIPPFloAccumulation(wippflo_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res,debug_cell)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(wippflo_auxvar_type) :: wippflo_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  PetscBool :: debug_cell
  
  PetscInt :: icomp, iphase
  
  PetscReal :: porosity
  PetscReal :: volume_over_dt
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use wippflo_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = wippflo_auxvar%effective_porosity
  
  ! accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] * 
    Res(iphase) = Res(iphase) + wippflo_auxvar%sat(iphase) * &
                                wippflo_auxvar%den(iphase)
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nflowspec) = Res(1:option%nflowspec) * &
                            porosity * volume_over_dt

end subroutine WIPPFloAccumulation

! ************************************************************************** !

subroutine WIPPFloFluxHarmonicPermOnly(wippflo_auxvar_up,global_auxvar_up, &
                                   material_auxvar_up, &
                                   wippflo_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   area, dist, upwind_direction_, &
                                   wippflo_parameter, &
                                   option,v_darcy,Res, &
                                   derivative_call, &
                                   update_upwind_direction_, &
                                   count_upwind_direction_flip_, &
                                   debug_connection)
  ! 
  ! Computes the internal flux terms for the residual based on harmonic 
  ! intrinsic permeability
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  use Upwind_Direction_module
  
  implicit none
  
  type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Res(option%nflowdof)
  PetscBool :: derivative_call
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscInt :: wat_comp_id, air_comp_id
  PetscInt :: iphase
  PetscBool :: upwind
  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, q
  PetscReal :: tot_mole_flux, wat_mole_flux, air_mole_flux
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  
  PetscReal :: temp_perm_up, temp_perm_dn

  PetscReal :: up_scale, dn_scale
  PetscReal :: dummy
  PetscInt :: iabs_upwind_direction1
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id

  Res = 0.d0
  v_darcy = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
!geh: we do not want to use the dot product with the unit vector, instead
!     use the principle direction stored in the upwind direction array
!  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
!  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  iabs_upwind_direction1 = iabs(upwind_direction_(1))
  select case(iabs_upwind_direction1)
    case(X_DIRECTION)
      perm_up = material_auxvar_up%permeability(perm_xx_index) 
      perm_dn = material_auxvar_dn%permeability(perm_xx_index) 
    case(Y_DIRECTION)
      perm_up = material_auxvar_up%permeability(perm_yy_index) 
      perm_dn = material_auxvar_dn%permeability(perm_yy_index) 
    case(Z_DIRECTION)
      perm_up = material_auxvar_up%permeability(perm_zz_index) 
      perm_dn = material_auxvar_dn%permeability(perm_zz_index) 
  end select
  
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_up%fracture) .and. &
      wippflo_use_fracture) then
    if (material_auxvar_up%fracture%vector(iabs_upwind_direction1) > 0.d0) then
      perm_up = perm_up * wippflo_auxvar_up%fracture_perm_scaling_factor
    endif
  endif
  if (associated(material_auxvar_dn%fracture) .and. &
      wippflo_use_fracture) then
    if (material_auxvar_dn%fracture%vector(iabs_upwind_direction1) > 0.d0) then
      perm_dn = perm_dn * wippflo_auxvar_dn%fracture_perm_scaling_factor
    endif
  endif
  
  perm_ave_over_dist(1) = (perm_up * perm_dn) / &
                          (dist_up*perm_dn + dist_dn*perm_up)
  temp_perm_up = wippflo_auxvar_up% &
                   klinkenberg_scaling_factor(iabs_upwind_direction1)*perm_up
  temp_perm_dn = wippflo_auxvar_dn% &
                   klinkenberg_scaling_factor(iabs_upwind_direction1)*perm_dn
  perm_ave_over_dist(2) = (temp_perm_up * temp_perm_dn) / &
                          (dist_up*temp_perm_dn + dist_dn*temp_perm_up)
      
  iphase = LIQUID_PHASE
  if (wippflo_auxvar_up%mobility(iphase) + &
      wippflo_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = WIPPFloAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           wippflo_auxvar_up%den_kg, &
                                           wippflo_auxvar_dn%den_kg, &
                                           dummy,dummy)
    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = wippflo_auxvar_up%pres(iphase) - &
                     wippflo_auxvar_dn%pres(iphase) + &
                     gravity_term
    up_scale = 0.d0
    dn_scale = 0.d0
    upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                             derivative_call, &
                             count_upwind_direction_flip_, &
                             liq_upwind_flip_count_by_res, &
                             liq_upwind_flip_count_by_jac)
    if (upwind) then
      up_scale = 1.d0
      mobility = wippflo_auxvar_up%mobility(iphase)
    else
      dn_scale = 1.d0
      mobility = wippflo_auxvar_dn%mobility(iphase)
    endif      

    if (mobility > floweps ) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = WIPPFloAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          wippflo_auxvar_up%den, &
                                          wippflo_auxvar_dn%den, &
                                          dummy,dummy)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      tot_mole_flux = q*density_ave
      ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
      !                                 xmol[kmol comp/kmol phase]
      wat_mole_flux = tot_mole_flux
      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      
    endif                   
  endif

  iphase = GAS_PHASE
  if (wippflo_auxvar_up%mobility(iphase) + &
      wippflo_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = WIPPFloAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           wippflo_auxvar_up%den_kg, &
                                           wippflo_auxvar_dn%den_kg, &
                                           dummy,dummy)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = wippflo_auxvar_up%pres(iphase) - &
                     wippflo_auxvar_dn%pres(iphase) + &
                     gravity_term
    ! if a gas phase does not exist on either side of the connection, the gas
    ! phase properties from the opposite side are used.
    up_scale = 0.d0
    dn_scale = 0.d0
    upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                             derivative_call, &
                             count_upwind_direction_flip_, &
                             gas_upwind_flip_count_by_res, &
                             gas_upwind_flip_count_by_jac)
    if (upwind) then
      up_scale = 1.d0
      mobility = wippflo_auxvar_up%mobility(iphase)
    else
      dn_scale = 1.d0
      mobility = wippflo_auxvar_dn%mobility(iphase)
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure
      density_ave = WIPPFloAverageDensity(iphase, &
                                          global_auxvar_up%istate, &
                                          global_auxvar_dn%istate, &
                                          wippflo_auxvar_up%den, &
                                          wippflo_auxvar_dn%den, &
                                          dummy,dummy)
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      tot_mole_flux = q*density_ave
      ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
      !                                 xmol[kmol comp/kmol phase]
      air_mole_flux = tot_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux

    endif               
  endif

end subroutine WIPPFloFluxHarmonicPermOnly

! ************************************************************************** !

subroutine WIPPFloFluxLumpedHarmonic(wippflo_auxvar_up,global_auxvar_up, &
                                     material_auxvar_up, &
                                     wippflo_auxvar_dn,global_auxvar_dn, &
                                     material_auxvar_dn, &
                                     area, dist, upwind_direction_, &
                                     wippflo_parameter, &
                                     option,v_darcy,Res, &
                                     derivative_call, &
                                     update_upwind_direction_, &
                                     count_upwind_direction_flip_, &
                                     debug_connection)
  ! 
  ! Computes the internal flux terms for the residual using a harmonic mean
  ! on the lumped density * permeability / viscosity.  The approach used in 
  ! BRAGFLO
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module
  use Material_Aux_class
  use Upwind_Direction_module
  
  implicit none
  
  type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area   ! area here is really area / alpha
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Res(option%nflowdof)
  PetscBool :: derivative_call
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_
  PetscBool :: debug_connection

  PetscReal :: upweight
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscInt :: wat_comp_id, air_comp_id
  PetscInt :: iphase
  PetscBool :: upwind
  
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: perm_rho_mu_area_ave_over_dist(2)
  PetscReal :: area_up, area_dn, area_ave
  PetscReal :: perm_up, perm_dn
  PetscReal :: dummy
  PetscReal :: delta_pressure
  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: rel_perm, q
  PetscReal :: wat_mole_flux, air_mole_flux
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  
  PetscReal :: temp_perm_up, temp_perm_dn

  PetscReal :: up_scale, dn_scale
  PetscInt :: iabs_upwind_direction1
  PetscReal :: perm_rho_mu_area_up(2), perm_rho_mu_area_dn(2)
  PetscReal :: gravity_
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id

  Res = 0.d0
  v_darcy = 0.d0
  ! to ensure that test values are set.
  perm_dn = 0.d0
  perm_dn = 1.d0/perm_dn
  perm_dn = 0.d0*perm_dn
  perm_up = perm_dn

  dist_up = dist(0)*dist(-1)
  dist_dn = dist(0)-dist_up
  upweight = 1.d0-dist(-1)
  gravity_ = option%gravity(3)
  dist_gravity = gravity_ * &
                 (wippflo_auxvar_dn%elevation - wippflo_auxvar_up%elevation)
!geh: we do not want to use the dot product with the unit vector, instead
!     use the principle direction stored in the upwind direction array
!  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
!  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  iabs_upwind_direction1 = iabs(upwind_direction_(1))
  select case(iabs_upwind_direction1)
    case(X_DIRECTION)
      perm_up = material_auxvar_up%permeability(perm_xx_index) 
      perm_dn = material_auxvar_dn%permeability(perm_xx_index) 
    case(Y_DIRECTION)
      perm_up = material_auxvar_up%permeability(perm_yy_index) 
      perm_dn = material_auxvar_dn%permeability(perm_yy_index) 
    case(Z_DIRECTION)
      perm_up = material_auxvar_up%permeability(perm_zz_index) 
      perm_dn = material_auxvar_dn%permeability(perm_zz_index) 
  end select
  
  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_up%fracture) .and. &
      wippflo_use_fracture) then
    if (material_auxvar_up%fracture%vector(iabs_upwind_direction1) > 0.d0) then
      perm_up = perm_up * wippflo_auxvar_up%fracture_perm_scaling_factor
    endif
  endif
  if (associated(material_auxvar_dn%fracture) .and. &
      wippflo_use_fracture) then
    if (material_auxvar_dn%fracture%vector(iabs_upwind_direction1) > 0.d0) then
      perm_dn = perm_dn * wippflo_auxvar_dn%fracture_perm_scaling_factor
    endif
  endif
  
  ! this scale by area
  area_up = wippflo_auxvar_up%alpha * area
  area_dn = wippflo_auxvar_dn%alpha * area
  area_ave = 0.5*(area_up+area_dn)
  perm_rho_mu_area_up(:) = perm_up * wippflo_auxvar_up%den / &
                           wippflo_auxvar_up%mu * area_up
  perm_rho_mu_area_dn(:) = perm_dn * wippflo_auxvar_dn%den / &
                           wippflo_auxvar_dn%mu * area_dn
  perm_rho_mu_area_up(2) = perm_rho_mu_area_up(2) * wippflo_auxvar_up% &
                       klinkenberg_scaling_factor(iabs_upwind_direction1)
  perm_rho_mu_area_dn(2) = perm_rho_mu_area_dn(2) * wippflo_auxvar_dn% &
                       klinkenberg_scaling_factor(iabs_upwind_direction1)

  ! this is an array(2)
  perm_rho_mu_area_ave_over_dist = &
    (perm_rho_mu_area_up * perm_rho_mu_area_dn) / &
    (dist_up*perm_rho_mu_area_dn + dist_dn*perm_rho_mu_area_up)
      
  iphase = LIQUID_PHASE
  density_kg_ave = upweight * wippflo_auxvar_up%den_kg(iphase) + &
                   (1.d0-upweight) * wippflo_auxvar_dn%den_kg(iphase)
  gravity_term = density_kg_ave * dist_gravity
  delta_pressure = wippflo_auxvar_up%pres(iphase) - &
                   wippflo_auxvar_dn%pres(iphase) + &
                   gravity_term

  up_scale = 0.d0
  dn_scale = 0.d0
  upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                           derivative_call, &
                           count_upwind_direction_flip_, &
                           liq_upwind_flip_count_by_res, &
                           liq_upwind_flip_count_by_jac)
  if (upwind) then
    up_scale = 1.d0
    rel_perm = wippflo_auxvar_up%kr(iphase)
  else
    dn_scale = 1.d0
    rel_perm = wippflo_auxvar_dn%kr(iphase)
  endif      

  ! wat_mole_flux[kmol/sec] = rho[kmol/m^3 phase] *
  !                           perm_area[m^4] / dist[m] * kr[-] / 
  !                           mu[Pa-sec] * dP[Pa]]
  wat_mole_flux = perm_rho_mu_area_ave_over_dist(iphase) * &
                  rel_perm * delta_pressure
  density_ave = 0.5d0*(wippflo_auxvar_up%den(iphase) + &
                       wippflo_auxvar_dn%den(iphase))
  ! v_darcy[m/sec] = wat_mole_flux[kmol/sec] / rho[kmol/m^3 phase] / 
  !                  area [m^2]
  v_darcy(iphase) = wat_mole_flux / density_ave / area_ave
  Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
  if (debug_connection) then
!    write(*,'("liq-t1X: ",9es12.4)') &
!      perm_rho_mu_area_ave_over_dist(iphase)*fmw_comp(iphase)* &
!      (dist_up+dist_dn)/area!, &
!      wippflo_auxvar_up%den(iphase)*fmw_comp(iphase), &
!      wippflo_auxvar_up%mu(iphase), perm_up, &
!      wippflo_auxvar_up%alpha, &
!      wippflo_auxvar_dn%den(iphase)*fmw_comp(iphase), &
!      wippflo_auxvar_dn%mu(iphase), perm_dn, &
!      wippflo_auxvar_dn%alpha
    write(*,'(9x,"a1l-up/dn: ",8es12.4)') &
      ! delta L = volume / (area*alpha)
      ! 1/(distL * deltaL * alpha)
      1.d0/(dist(0)*material_auxvar_up%volume/area), &
      1.d0/(dist(0)*material_auxvar_dn%volume/area)
    write(*,'(9x,"liq-up,t1l,dp,t2l: ",l1,x,8es12.4)') upwind, &
      perm_rho_mu_area_ave_over_dist(iphase)*fmw_comp(iphase)* &
      (dist_up+dist_dn)/area, &
      delta_pressure, &
      rel_perm!, &
!      gravity_term,
!      wippflo_auxvar_up%pres(iphase), &
!      wippflo_auxvar_dn%pres(iphase)
  endif

  iphase = GAS_PHASE
  density_kg_ave = upweight * wippflo_auxvar_up%den_kg(iphase) + &
                   (1.d0-upweight) * wippflo_auxvar_dn%den_kg(iphase)
  gravity_term = density_kg_ave * dist_gravity
  delta_pressure = wippflo_auxvar_up%pres(iphase) - &
                   wippflo_auxvar_dn%pres(iphase) + &
                   gravity_term

  ! if a gas phase does not exist on either side of the connection, the gas
  ! phase properties from the opposite side are used.
  up_scale = 0.d0
  dn_scale = 0.d0
  upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                           derivative_call, &
                           count_upwind_direction_flip_, &
                           gas_upwind_flip_count_by_res, &
                           gas_upwind_flip_count_by_jac)
  if (upwind) then
    up_scale = 1.d0
    rel_perm = wippflo_auxvar_up%kr(iphase)
  else
    dn_scale = 1.d0
    rel_perm = wippflo_auxvar_dn%kr(iphase)
  endif      

  ! air_mole_flux[kmol/sec] = rho[kmol/m^3 phase] *
  !                           perm_area[m^4] / dist[m] * kr[-] / 
  !                           mu[Pa-sec] * dP[Pa]]
  air_mole_flux = perm_rho_mu_area_ave_over_dist(iphase) * &
                  rel_perm * delta_pressure
  density_ave = 0.5d0*(wippflo_auxvar_up%den(iphase) + &
                       wippflo_auxvar_dn%den(iphase))
  ! v_darcy[m/sec] = air_mole_flux[kmol/sec] / rho[kmol/m^3 phase] / 
  !                  area [m^2]
  v_darcy(iphase) = air_mole_flux / density_ave / area_ave
  Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
  if (debug_connection) then
    write(*,'(9x,"gas-up,t1l,dp,t2l: ",l1,x,8es12.4)') upwind, &
      perm_rho_mu_area_ave_over_dist(iphase)*fmw_comp(iphase)* &
      (dist_up+dist_dn)/area, &
      delta_pressure, &
      rel_perm!, &
!      gravity_term,
!      wippflo_auxvar_up%pres(iphase), &
!      wippflo_auxvar_dn%pres(iphase)
  endif

end subroutine WIPPFloFluxLumpedHarmonic

! ************************************************************************** !

subroutine WIPPFloBCFluxHarmonicPermOnly(ibndtype,auxvar_mapping,auxvars, &
                                     wippflo_auxvar_up,global_auxvar_up, &
                                     wippflo_auxvar_dn,global_auxvar_dn, &
                                     material_auxvar_dn, &
                                     area,dist,upwind_direction_, &
                                     wippflo_parameter, &
                                     option,v_darcy,Res, &
                                     derivative_call, &
                                     update_upwind_direction_, &
                                     count_upwind_direction_flip_, &
                                     debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module                              
  use Material_Aux_class
  use Upwind_Direction_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(WIPPFLO_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: Res(1:option%nflowdof)
  PetscBool :: derivative_call
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_
  PetscBool :: debug_connection
  
  PetscInt :: wat_comp_id, air_comp_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: mobility, q 
  PetscReal :: tot_mole_flux
  PetscReal :: perm_dn
  PetscReal :: boundary_pressure
  PetscReal :: tempreal
  PetscReal :: wat_mole_flux, air_mole_flux
  PetscBool :: upwind
  PetscInt :: iabs_upwind_direction1
  PetscReal :: dummy
  PetscReal :: dn_scale

  PetscInt :: idof
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id

  Res = 0.d0
  v_darcy = 0.d0  

!geh: we do not want to use the dot product with the unit vector, instead
!     use the principle direction stored in the upwind direction array
!  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  iabs_upwind_direction1 = iabs(upwind_direction_(1))
  select case(iabs_upwind_direction1)
    case(X_DIRECTION)
      perm_dn = material_auxvar_dn%permeability(perm_xx_index) 
    case(Y_DIRECTION)
      perm_dn = material_auxvar_dn%permeability(perm_yy_index) 
    case(Z_DIRECTION)
      perm_dn = material_auxvar_dn%permeability(perm_zz_index) 
  end select

  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_dn%fracture) .and. &
      wippflo_use_fracture) then
    if (material_auxvar_dn%fracture%vector(iabs_upwind_direction1) > 0.d0) then
      perm_dn = perm_dn * wippflo_auxvar_dn%fracture_perm_scaling_factor
    endif
  endif
  
  perm_dn_adj(1) = perm_dn
  perm_dn_adj(2) = wippflo_auxvar_dn% &
                     klinkenberg_scaling_factor(iabs_upwind_direction1)*perm_dn
  
  iphase = LIQUID_PHASE
  mobility = 0.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      if (wippflo_auxvar_up%mobility(iphase) + &
          wippflo_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(WIPPFLO_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(WIPPFLO_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        boundary_pressure = wippflo_auxvar_up%pres(iphase)
        density_kg_ave = WIPPFloAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                wippflo_auxvar_up%den_kg, &
                                                wippflo_auxvar_dn%den_kg, &
                                                dummy,dummy)
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          wippflo_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              wippflo_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
          endif
        endif
        dn_scale = 0.d0
        upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                                 derivative_call, &
                                 count_upwind_direction_flip_, &
                                 liq_bc_upwind_flip_count_by_res, &
                                 liq_bc_upwind_flip_count_by_jac)
        if (upwind) then
          mobility = wippflo_auxvar_up%mobility(iphase)
        else
          dn_scale = 1.d0        
          mobility = wippflo_auxvar_dn%mobility(iphase)
        endif      

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = WIPPFloAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            wippflo_auxvar_up%den, &
                                            wippflo_auxvar_dn%den, &
                                            dummy,dummy)
      endif
    case(NEUMANN_BC)
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(WIPPFLO_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(WIPPFLO_GAS_FLUX_INDEX)
      end select
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = wippflo_auxvar_up%den(iphase)
        else 
          dn_scale = 1.d0
          density_ave = wippflo_auxvar_dn%den(iphase)
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in WIPPFloBCFlux phase loop.'
      call PrintErrMsg(option)
  end select
  if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    wat_mole_flux = tot_mole_flux
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
   
  endif                   

  iphase = GAS_PHASE
  mobility = 0.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,&
         HYDROSTATIC_CONDUCTANCE_BC)
      if (wippflo_auxvar_up%mobility(iphase) + &
          wippflo_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(WIPPFLO_LIQUID_CONDUCTANCE_INDEX)
            case(GAS_PHASE)
              idof = auxvar_mapping(WIPPFLO_GAS_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        boundary_pressure = wippflo_auxvar_up%pres(iphase)
        density_kg_ave = WIPPFloAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                wippflo_auxvar_up%den_kg, &
                                                wippflo_auxvar_dn%den_kg, &
                                                dummy,dummy)
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          wippflo_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              wippflo_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
          endif
        endif
        dn_scale = 0.d0
        ! don't expect the derivative to match precisely at delta_pressure = 0
        ! due to potential switch in direction for numerically perturbed
        ! residual
        upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                                 derivative_call, &
                                 count_upwind_direction_flip_, &
                                 gas_bc_upwind_flip_count_by_res, &
                                 gas_bc_upwind_flip_count_by_jac)
        if (upwind) then
          mobility = wippflo_auxvar_up%mobility(iphase)
        else
          dn_scale = 1.d0        
          mobility = wippflo_auxvar_dn%mobility(iphase)
        endif      
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = WIPPFloAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            wippflo_auxvar_up%den, &
                                            wippflo_auxvar_dn%den, &
                                            dummy,dummy)    
      endif
    case(NEUMANN_BC)
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(WIPPFLO_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(WIPPFLO_GAS_FLUX_INDEX)
      end select
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = wippflo_auxvar_up%den(iphase)
        else 
          dn_scale = 1.d0
          density_ave = wippflo_auxvar_dn%den(iphase)
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in WIPPFloBCFlux phase loop.'
      call PrintErrMsg(option)
  end select

  if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    air_mole_flux = tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      
  endif                   

end subroutine WIPPFloBCFluxHarmonicPermOnly

! ************************************************************************** !

subroutine WIPPFloBCFluxLumpedHarmonic(ibndtype,auxvar_mapping,auxvars, &
                                       wippflo_auxvar_up,global_auxvar_up, &
                                       wippflo_auxvar_dn,global_auxvar_dn, &
                                       material_auxvar_dn, &
                                       area,dist,upwind_direction_, &
                                       wippflo_parameter, &
                                       option,v_darcy,Res, &
                                       derivative_call, &
                                       update_upwind_direction_, &
                                       count_upwind_direction_flip_, &
                                       debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module                              
  use Material_Aux_class
  use Upwind_Direction_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(WIPPFLO_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area ! this is the actual area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: Res(1:option%nflowdof)
  PetscBool :: derivative_call
  PetscBool :: update_upwind_direction_
  PetscBool :: count_upwind_direction_flip_
  PetscBool :: debug_connection
  
  PetscInt :: wat_comp_id, air_comp_id
  PetscInt :: icomp, iphase
  PetscInt :: bc_type
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: perm_dn_adj(2)
  PetscReal :: perm_ave_over_dist
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: dist_gravity
  PetscReal :: rel_perm, viscosity, q 
  PetscReal :: tot_mole_flux
  PetscReal :: perm_dn
  PetscReal :: boundary_pressure
  PetscReal :: tempreal
  PetscReal :: wat_mole_flux, air_mole_flux
  PetscBool :: upwind
  PetscInt :: iabs_upwind_direction1

  PetscReal :: dn_scale
  PetscReal :: dummy

  PetscInt :: idof
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id

  Res = 0.d0
  v_darcy = 0.d0  
  ! to ensure that test values are set.
  perm_dn = 0.d0
  perm_dn = 1.d0/perm_dn
  perm_dn = 0.d0*perm_dn
  density_ave = perm_dn

!geh: we do not want to use the dot product with the unit vector, instead
!     use the principle direction stored in the upwind direction array
!  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  iabs_upwind_direction1 = iabs(upwind_direction_(1))
  select case(iabs_upwind_direction1)
    case(X_DIRECTION)
      perm_dn = material_auxvar_dn%permeability(perm_xx_index) 
    case(Y_DIRECTION)
      perm_dn = material_auxvar_dn%permeability(perm_yy_index) 
    case(Z_DIRECTION)
      perm_dn = material_auxvar_dn%permeability(perm_zz_index) 
  end select

  ! Fracture permeability change only available for structured grid (Heeho)
  if (associated(material_auxvar_dn%fracture) .and. &
      wippflo_use_fracture) then
    if (material_auxvar_dn%fracture%vector(iabs_upwind_direction1) > 0.d0) then
      perm_dn = perm_dn * wippflo_auxvar_dn%fracture_perm_scaling_factor
    endif
  endif
  
  perm_dn_adj(1) = perm_dn
  perm_dn_adj(2) = wippflo_auxvar_dn% &
                     klinkenberg_scaling_factor(iabs_upwind_direction1)*perm_dn
  
  iphase = LIQUID_PHASE
  rel_perm = 0.d0
  bc_type = ibndtype(iphase)

  ! dist(0) = scalar - magnitude of distance
  ! gravity = vector(3)
  ! dist(1:3) = vector(3) - unit vector
  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

  select case(bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
        select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(WIPPFLO_LIQUID_CONDUCTANCE_INDEX)
          case(GAS_PHASE)
            idof = auxvar_mapping(WIPPFLO_GAS_CONDUCTANCE_INDEX)
        end select        
        perm_ave_over_dist = auxvars(idof)
      else
        perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
      endif
      
      boundary_pressure = wippflo_auxvar_up%pres(iphase)
      !geh: use density and viscosity at boundary
      gravity_term = wippflo_auxvar_up%den_kg(iphase) * dist_gravity
      viscosity = wippflo_auxvar_up%mu(iphase)
      delta_pressure = boundary_pressure - &
                       wippflo_auxvar_dn%pres(iphase) + &
                       gravity_term
      if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
          bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
            ! flow in         ! boundary cell is <= pref
        if (delta_pressure > 0.d0 .and. &
            wippflo_auxvar_up%pres(iphase) - &
              option%reference_pressure < eps) then
          delta_pressure = 0.d0
        endif
      endif
      dn_scale = 0.d0
      upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                               derivative_call, &
                               count_upwind_direction_flip_, &
                               liq_bc_upwind_flip_count_by_res, &
                               liq_bc_upwind_flip_count_by_jac)
      if (upwind) then
        rel_perm = wippflo_auxvar_up%kr(iphase)
      else
        dn_scale = 1.d0        
        rel_perm = wippflo_auxvar_dn%kr(iphase)
      endif      

      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist * rel_perm / viscosity * &
                        delta_pressure
      ! density at interface
      density_ave = wippflo_auxvar_up%den(iphase)

    case(NEUMANN_BC)
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(WIPPFLO_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(WIPPFLO_GAS_FLUX_INDEX)
      end select
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = wippflo_auxvar_up%den(iphase)
        else 
          dn_scale = 1.d0
          density_ave = wippflo_auxvar_dn%den(iphase)
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in BRAGFloBCFlux phase loop.'
      call PrintErrMsg(option)
  end select
  if (dabs(v_darcy(iphase)) > 0.d0 .or. rel_perm > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    wat_mole_flux = tot_mole_flux
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
  endif                   

  iphase = GAS_PHASE
  rel_perm = 0.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
        select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(WIPPFLO_LIQUID_CONDUCTANCE_INDEX)
          case(GAS_PHASE)
            idof = auxvar_mapping(WIPPFLO_GAS_CONDUCTANCE_INDEX)
        end select        
        perm_ave_over_dist = auxvars(idof)
      else
        perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
      endif
      
      boundary_pressure = wippflo_auxvar_up%pres(iphase)
      !geh: use density and viscosity at boundary
      gravity_term = wippflo_auxvar_up%den_kg(iphase) * dist_gravity
      viscosity = wippflo_auxvar_up%mu(iphase)
      delta_pressure = boundary_pressure - &
                       wippflo_auxvar_dn%pres(iphase) + &
                       gravity_term
      if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
          bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
            ! flow in         ! boundary cell is <= pref
        if (delta_pressure > 0.d0 .and. &
            wippflo_auxvar_up%pres(iphase) - &
              option%reference_pressure < eps) then
          delta_pressure = 0.d0
        endif
      endif
      dn_scale = 0.d0
      ! don't expect the derivative to match precisely at delta_pressure = 0
      ! due to potential switch in direction for numerically perturbed
      ! residual
      upwind = UpwindDirection(upwind_direction_(iphase),delta_pressure, &
                               derivative_call, &
                               count_upwind_direction_flip_, &
                               gas_bc_upwind_flip_count_by_res, &
                               gas_bc_upwind_flip_count_by_jac)
      if (upwind) then
        rel_perm = wippflo_auxvar_up%kr(iphase)
      else
        dn_scale = 1.d0        
        rel_perm = wippflo_auxvar_dn%kr(iphase)
      endif      
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist * rel_perm / viscosity * &
                        delta_pressure
      ! density at interface
      density_ave = wippflo_auxvar_up%den(iphase)

    case(NEUMANN_BC)
      dn_scale = 0.d0
      select case(iphase)
        case(LIQUID_PHASE)
          idof = auxvar_mapping(WIPPFLO_LIQUID_FLUX_INDEX)
        case(GAS_PHASE)
          idof = auxvar_mapping(WIPPFLO_GAS_FLUX_INDEX)
      end select
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(iphase) = auxvars(idof)
        if (v_darcy(iphase) > 0.d0) then 
          density_ave = wippflo_auxvar_up%den(iphase)
        else 
          dn_scale = 1.d0
          density_ave = wippflo_auxvar_dn%den(iphase)
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in BRAGFloBCFlux phase loop.'
      call PrintErrMsg(option)
  end select

  if (dabs(v_darcy(iphase)) > 0.d0 .or. rel_perm > 0.d0) then
    ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(iphase) * area  
    ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
    !                             density_ave[kmol phase/m^3 phase]        
    tot_mole_flux = q*density_ave
    ! comp_mole_flux[kmol comp/sec] = tot_mole_flux[kmol phase/sec] * 
    !                                 xmol[kmol comp/kmol phase]
    air_mole_flux = tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
  endif                   

end subroutine WIPPFloBCFluxLumpedHarmonic

! ************************************************************************** !

subroutine WIPPFloSrcSink(option,qsrc,flow_src_sink_type, &
                          wippflo_auxvar,global_auxvar,material_auxvar, &
                          ss_flow_vol_flux,scale,Res,debug_cell)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module
  use Material_Aux_class
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(wippflo_auxvar_type) :: wippflo_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscBool :: debug_cell
      
  PetscReal :: qsrc_mol
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscReal :: dden_bool
  PetscErrorCode :: ierr

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  
  Res = 0.d0
  
  ! liquid phase
  qsrc_mol = 0.d0
  dden_bool = 0.d0
  select case(flow_src_sink_type)
    case(MASS_RATE_SS)
      qsrc_mol = qsrc(wat_comp_id)/fmw_comp(wat_comp_id) ! kg/sec -> kmol/sec
    case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
      qsrc_mol = qsrc(wat_comp_id)/fmw_comp(wat_comp_id)*scale 
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_mol = qsrc(wat_comp_id)*wippflo_auxvar%den(wat_comp_id) ! den = kmol/m^3
      dden_bool = 1.d0
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_mol = qsrc(wat_comp_id)*wippflo_auxvar%den(wat_comp_id)*scale
      dden_bool = 1.d0
  end select
  ss_flow_vol_flux(wat_comp_id) = qsrc_mol/wippflo_auxvar%den(wat_comp_id)
  Res(wat_comp_id) = qsrc_mol

  ! gas phase
  qsrc_mol = 0.d0
  dden_bool = 0.d0
  select case(flow_src_sink_type)
    case(MASS_RATE_SS)
      qsrc_mol = qsrc(air_comp_id)/fmw_comp(air_comp_id) ! kg/sec -> kmol/sec
    case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
      qsrc_mol = qsrc(air_comp_id)/fmw_comp(air_comp_id)*scale 
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_mol = qsrc(air_comp_id)*wippflo_auxvar%den(air_comp_id) ! den = kmol/m^3
      dden_bool = 1.d0
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_mol = qsrc(air_comp_id)*wippflo_auxvar%den(air_comp_id)*scale
      dden_bool = 1.d0
  end select
  ss_flow_vol_flux(air_comp_id) = qsrc_mol/wippflo_auxvar%den(air_comp_id)
  Res(air_comp_id) = qsrc_mol

  if (dabs(qsrc(TWO_INTEGER)) < 1.d-40 .and. &
      qsrc(ONE_INTEGER) < 0.d0) then ! extraction only
    Res(TWO_INTEGER) = qsrc_mol
    ss_flow_vol_flux(air_comp_id) = qsrc_mol/wippflo_auxvar%den(TWO_INTEGER)
  endif

end subroutine WIPPFloSrcSink

! ************************************************************************** !

subroutine WIPPFloAccumDerivative(wippflo_auxvar,global_auxvar, &
                                  material_auxvar, soil_heat_capacity, &
                                  option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none

  type(wippflo_auxvar_type) :: wippflo_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  call WIPPFloAccumulation(wippflo_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar,soil_heat_capacity,option, &
                           res,PETSC_FALSE)
                           
  do idof = 1, option%nflowdof
    call WIPPFloAccumulation(wippflo_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar,soil_heat_capacity, &
                             option,res_pert,PETSC_FALSE)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvar(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine WIPPFloAccumDerivative

! ************************************************************************** !

subroutine XXFluxDerivative(wippflo_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 wippflo_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 area, dist, &
                                 upwind_direction_, &
                                 wippflo_parameter, &
                                 option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module
  use Material_Aux_class
  use Upwind_Direction_module, only : count_upwind_direction_flip
  
  implicit none
  
  type(wippflo_auxvar_type) :: wippflo_auxvar_up(0:), wippflo_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res_up(option%nflowdof), res_dn(option%nflowdof)
  PetscReal :: res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0

  option%iflag = -2
  call XXFlux(wippflo_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   wippflo_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   area,dist,upwind_direction_, &
                   wippflo_parameter, &
                   option,v_darcy,res_up, &
                   PETSC_TRUE, & ! derivative call 
                   PETSC_FALSE, & ! update the upwind direction
                   ! avoid double counting upwind direction flip
                   PETSC_FALSE, & ! count upwind direction flip
                   PETSC_FALSE)
  res_dn = res_up

  if (wippflo_jacobian_test) then
    print *, 'res_dn: ', res_dn
    print *, 'res_up: ', res_up
  endif
 
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call XXFlux(wippflo_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up, &
                     wippflo_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist,upwind_direction_, &
                     wippflo_parameter, &
                     option,v_darcy,res_pert, &
                     PETSC_TRUE, & ! derivative call
                     PETSC_FALSE, & ! update the upwind direction
                     count_upwind_direction_flip, & 
                     PETSC_FALSE)
    if (wippflo_jacobian_test) then
      if (wippflo_jacobian_test_xdof > 0 .and. &
          idof == 2-mod(wippflo_jacobian_test_xdof,2)) then
        print *, 'res_pert_up: ', res_pert
      endif
    endif
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res_up(irow)) / &
                       wippflo_auxvar_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call XXFlux(wippflo_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                     material_auxvar_up, &
                     wippflo_auxvar_dn(idof),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist,upwind_direction_, &
                     wippflo_parameter, &
                     option,v_darcy,res_pert, &
                     PETSC_TRUE, & ! derivative call
                     PETSC_FALSE, & ! update the upwind direction
                     count_upwind_direction_flip, & 
                     PETSC_FALSE)
    if (wippflo_jacobian_test) then
      if (wippflo_jacobian_test_xdof > 0 .and. &
          idof == 2-mod(wippflo_jacobian_test_xdof,2)) then
        print *, 'res_pert_dn: ', res_pert
      endif
    endif
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res_dn(irow)) / &
                       wippflo_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine XXFluxDerivative

! ************************************************************************** !

subroutine XXBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   wippflo_auxvar_up, &
                                   global_auxvar_up, &
                                   wippflo_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   area,dist,upwind_direction_, &
                                   wippflo_parameter, &
                                   option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module 
  use Material_Aux_class
  use Upwind_Direction_module, only : count_upwind_direction_flip
  
  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(WIPPFLO_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction_(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jdn = 0.d0

  option%iflag = -2
  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     wippflo_auxvar_up,global_auxvar_up, &
                     wippflo_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist,upwind_direction_, &
                     wippflo_parameter, &
                     option,v_darcy,res, &
                     PETSC_TRUE, & ! derivative call
                     PETSC_FALSE, & ! update the upwind direction
                     ! avoid double counting upwind direction flip
                     PETSC_FALSE, & ! count upwind direction flip
                     PETSC_FALSE)

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       wippflo_auxvar_up,global_auxvar_up, &
                       wippflo_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,upwind_direction_, &
                       wippflo_parameter, &
                       option,v_darcy,res_pert, &
                       PETSC_TRUE, & ! derivative call
                       PETSC_FALSE, & ! update the upwind direction
                       count_upwind_direction_flip, & 
                       PETSC_FALSE)   
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

end subroutine XXBCFluxDerivative

! ************************************************************************** !

subroutine WIPPFloSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                    wippflo_auxvars,global_auxvar, &
                                    material_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(wippflo_auxvar_type) :: wippflo_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3
  ! unperturbed wippflo_auxvars value
  call WIPPFloSrcSink(option,qsrc,flow_src_sink_type, &
                      wippflo_auxvars(ZERO_INTEGER),global_auxvar, &
                      material_auxvar,dummy_real,scale,res,PETSC_FALSE)
                      
  ! perturbed wippflo_auxvars values
  do idof = 1, option%nflowdof
    call WIPPFloSrcSink(option,qsrc,flow_src_sink_type, &
                        wippflo_auxvars(idof),global_auxvar, &
                        material_auxvar,dummy_real, &
                        scale,res_pert,PETSC_FALSE)            
    do irow = 1, option%nflowdof
      Jac(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvars(idof)%pert
    enddo !irow
  enddo ! idof
  
end subroutine WIPPFloSrcSinkDerivative

! ************************************************************************** !

function WIPPFloAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: WIPPFloAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  WIPPFloAverageDensity = 0.d0
  if (iphase == LIQUID_PHASE) then
    WIPPFloAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
    dden_up = 0.5d0
    dden_dn = 0.5d0
  else if (iphase == GAS_PHASE) then
    WIPPFloAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
    dden_up = 0.5d0
    dden_dn = 0.5d0      
  endif

end function WIPPFloAverageDensity

end module WIPP_Flow_Common_module
