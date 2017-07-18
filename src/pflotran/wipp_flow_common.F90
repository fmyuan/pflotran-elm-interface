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

  ! variables that track the number of times the upwind direction changes
  ! during the residual and Jacobian calculations.
  PetscInt, public :: liq_upwind_flip_count_by_res
  PetscInt, public :: gas_upwind_flip_count_by_res
  PetscInt, public :: liq_bc_upwind_flip_count_by_res
  PetscInt, public :: gas_bc_upwind_flip_count_by_res
  PetscInt, public :: liq_upwind_flip_count_by_jac
  PetscInt, public :: gas_upwind_flip_count_by_jac
  PetscInt, public :: liq_bc_upwind_flip_count_by_jac
  PetscInt, public :: gas_bc_upwind_flip_count_by_jac

  public :: WIPPFloAccumulation, &
            WIPPFloFlux, &
            WIPPFloBCFlux, &
            WIPPFloSrcSink, &
            WIPPFloAccumDerivative, &
            WIPPFloFluxDerivative, &
            WIPPFloBCFluxDerivative, &
            WIPPFloSrcSinkDerivative
            
  public :: WIPPFloDiffJacobian, &
            WIPPFloAuxVarDiff

contains

! ************************************************************************** !

subroutine WIPPFloAccumulation(wippflo_auxvar,global_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res,Jac, &
                               analytical_derivatives,debug_cell)
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
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscBool :: analytical_derivatives
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
  
  if (analytical_derivatives) then
    Jac = 0.d0
    ! ----------
    ! Water Equation
    ! por * (satl * denl * Xwl + satg * deng * Xwg)
    ! ---
    ! w/respect to gas pressure
    ! dpor_dp * (satl * denl * Xwl + satg * deng * Xwg) +
    ! por * (satl * ddenl_dpg * Xwl + satg * ddeng_dpg * Xwg) + 
    ! por * (satl * denl * dXwl_dpg + satg * deng * dXwg_dpg)
    Jac(1,1) = &
      wippflo_auxvar%d%por_p * &
        (wippflo_auxvar%sat(1) * wippflo_auxvar%den(1)) + &
      porosity * &
                              ! denl_pl = denl_pg
        (wippflo_auxvar%sat(1) * wippflo_auxvar%d%denl_pl)
    Jac(1,2) = porosity * &
      (-1.d0 * wippflo_auxvar%den(1))
    ! ----------
    ! Gas Equation
    ! por * (satl * denl * Xal + satg * deng * Xag)
    ! ---
    ! w/respect to gas pressure
    ! dpor_dp * (satl * denl * Xal + satg * deng * Xag) +
    ! por * (satl * ddenl_dpg * Xal + satg * ddeng_dpg * Xag) + 
    ! por * (satl * denl * dXal_dpg + satg * deng * dXag_dpg)
    Jac(2,1) = &
      wippflo_auxvar%d%por_p * &
        (wippflo_auxvar%sat(2) * wippflo_auxvar%den(2)) + &
      porosity * &
        (wippflo_auxvar%sat(2) * wippflo_auxvar%d%deng_pg)
    ! w/respect to gas saturation
    ! porosity, density, mole fraction are independent of gas saturation
    ! por * (dsatl_dsatg * denl * Xal + dsatg_dsatg * deng * Xag)
    ! dsatl_dsatg = -1.
    ! dsatg_dsatg = 1.
    Jac(2,2) = porosity * &
      (1.d0 * wippflo_auxvar%den(2))
    Jac = Jac * volume_over_dt
  endif
  
end subroutine WIPPFloAccumulation

! ************************************************************************** !

subroutine WIPPFloFlux(wippflo_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       wippflo_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       area, dist, upwind_direction, &
                       wippflo_parameter, &
                       option,v_darcy,Res,Jup,Jdn, &
                       derivative_call, &
                       fix_upwind_direction, &
                       update_upwind_direction, &
                       count_upwind_direction_flip, &
                       analytical_derivatives, &
                       debug_connection)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  
  implicit none
  
  type(wippflo_auxvar_type) :: wippflo_auxvar_up, wippflo_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Res(option%nflowdof)
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscBool :: derivative_call
  PetscBool :: fix_upwind_direction
  PetscBool :: update_upwind_direction
  PetscBool :: count_upwind_direction_flip
  PetscBool :: analytical_derivatives
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

  ! Darcy flux
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dpaup, ddelta_pressure_dpadn
  
  PetscReal :: up_scale, dn_scale
  PetscReal :: tot_mole_flux_ddel_pressure
  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dtot_mole_flux_dp, dtot_mole_flux_dsatg
  PetscReal :: dpl_dsatg
  PetscReal :: ddelta_pressure_pl
  PetscInt :: prev_upwind_direction
  PetscInt :: new_upwind_direction
  PetscInt :: iabs_upwind_direction1
  
  ! DELETE
  
  PetscReal :: Jlup(2,2), Jldn(2,2)
  PetscReal :: Jgup(2,2), Jgdn(2,2)

  wat_comp_id = option%water_id
  air_comp_id = option%air_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
!geh: we do not want to use the dot product with the unit vector, instead
!     use the principle direction stored in the upwind direction array
!  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
!  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  iabs_upwind_direction1 = iabs(upwind_direction(1))
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
      
  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0  
  
  v_darcy = 0.d0

  !TODO(geh): merge phases and use arrays for partial derivatives.  This can 
  !           only be done after analytical derivatives are set up.
  iphase = LIQUID_PHASE
  if (wippflo_auxvar_up%mobility(iphase) + &
      wippflo_auxvar_dn%mobility(iphase) > eps) then
    
    density_kg_ave = WIPPFloAverageDensity(iphase, &
                                           global_auxvar_up%istate, &
                                           global_auxvar_dn%istate, &
                                           wippflo_auxvar_up%den_kg, &
                                           wippflo_auxvar_dn%den_kg, &
                                           ddensity_kg_ave_dden_kg_up, &
                                           ddensity_kg_ave_dden_kg_dn)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = wippflo_auxvar_up%pres(iphase) - &
                     wippflo_auxvar_dn%pres(iphase) + &
                     gravity_term
    if (analytical_derivatives) then
      ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_up * &
                             wippflo_auxvar_up%d%denl_pl * fmw_comp(iphase)
      ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_dn * &
                             wippflo_auxvar_dn%d%denl_pl * fmw_comp(iphase)
    endif
    up_scale = 0.d0
    dn_scale = 0.d0
    if (fix_upwind_direction) then
      if (update_upwind_direction .or. count_upwind_direction_flip) then
        prev_upwind_direction = upwind_direction(iphase)
        if (delta_pressure >= 0.d0) then
          ! positive means upstream
          new_upwind_direction = iabs(prev_upwind_direction)
        else
          ! negative means downstream
          new_upwind_direction = -iabs(prev_upwind_direction)
        endif 
        if (count_upwind_direction_flip) then
          if (new_upwind_direction /= prev_upwind_direction) then
            if (derivative_call) then
              liq_upwind_flip_count_by_jac = liq_upwind_flip_count_by_jac + 1
            else
              liq_upwind_flip_count_by_res = liq_upwind_flip_count_by_res + 1
            endif
          endif
        endif
        if (update_upwind_direction) then
          upwind_direction(iphase) = new_upwind_direction
        endif
      endif
      upwind = (upwind_direction(iphase) > 0)
    else
      upwind = (delta_pressure >= 0.d0)
    endif
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
      wat_mole_flux = tot_mole_flux
      Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
      
      if (analytical_derivatives) then
        option%io_buffer = 'Derivatives must be posed in terms of liquid &
          &pressure instead of gas pressure.'
        call printErrMsg(option)
        Jlup = 0.d0
        Jldn = 0.d0
        ! Upstream Cell
        ! derivative wrt gas pressure
        ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
        !   liquid pressure derivatives.
        ! derivative total mole flux wrt gas pressure
        dtot_mole_flux_dp = &
          ! ave. liquid density
          q * ddensity_ave_dden_up * wippflo_auxvar_up%d%denl_pl + &
          ! liquid mobility
          up_scale * &
          tot_mole_flux / mobility * wippflo_auxvar_up%d%mobilityl_pl + &
          ! pressure gradient
          tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
        ! derivative water wrt gas pressure
        Jlup(1,1) = dtot_mole_flux_dp
        ! derivative air wrt gas pressure
            
        ! derivative wrt gas saturation
        ! pl = pg - pc(satg)
        dpl_dsatg = -1.d0 * wippflo_auxvar_up%d%pc_satg
        ! delta pressure = plup - pldn
        ddelta_pressure_pl = 1.d0          
        ! derivative total mole flux wrt gas saturation
        dtot_mole_flux_dsatg = &
          ! liquid viscosity
          ! since liquid viscosity in a two phase state is a function
          ! of total pressure (gas pressure), there is not derivative
          ! wrt gas saturation
          !up_scale * &
          !tot_mole_flux / mobility * &
          !wippflo_auxvar_up%d%mobilityl_pl * dpl_dsatg + &
          ! relative permeability
          up_scale * &
          tot_mole_flux / mobility * &
          wippflo_auxvar_up%d%mobilityl_satg + &
          ! pressure gradient
          tot_mole_flux_ddel_pressure * ddelta_pressure_pl * dpl_dsatg         
        ! derivative water wrt gas saturation
        Jlup(1,2) = dtot_mole_flux_dsatg
        ! derivative air wrt gas saturation

        ! Downstream Cell
        ! derivative wrt gas pressure
        ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
        !   liquid pressure derivatives.
        ! derivative total mole flux wrt gas pressure
        dtot_mole_flux_dp = &
          ! ave. liquid density
          q * ddensity_ave_dden_dn *wippflo_auxvar_dn%d%denl_pl + &
          ! liquid mobility
          dn_scale * &
          tot_mole_flux / mobility * wippflo_auxvar_dn%d%mobilityl_pl + &
          ! pressure gradient
          tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
        ! derivative water wrt gas pressure
        Jldn(1,1) = dtot_mole_flux_dp
        ! derivative air wrt gas pressure
            
        ! derivative wrt gas saturation
        ! pl = pg - pc(satg)
        dpl_dsatg = -1.d0 * wippflo_auxvar_dn%d%pc_satg
        ! delta pressure = plup - pldn
        ddelta_pressure_pl = -1.d0
        ! derivative total mole flux wrt gas saturation
        dtot_mole_flux_dsatg = &
          ! liquid viscosity
          ! since liquid viscosity in a two phase state is a function
          ! of total pressure (gas pressure), there is not derivative
          ! wrt gas saturation
          !dn_scale * &
          !tot_mole_flux / mobility * &
          !wippflo_auxvar_dn%d%mobilityl_pl * dpl_dsatg + &
          ! relative permeability
          dn_scale * &
          tot_mole_flux / mobility * &
          wippflo_auxvar_dn%d%mobilityl_satg + &
          !pressure gradient
          tot_mole_flux_ddel_pressure * ddelta_pressure_pl * dpl_dsatg
        ! derivative water wrt gas saturation
        Jldn(1,2) = dtot_mole_flux_dsatg
        ! derivative air wrt gas saturation
          
        Jup = Jup + Jlup
        Jdn = Jdn + Jldn
      endif
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
                                           ddensity_kg_ave_dden_kg_up, &
                                           ddensity_kg_ave_dden_kg_dn)

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = wippflo_auxvar_up%pres(iphase) - &
                     wippflo_auxvar_dn%pres(iphase) + &
                     gravity_term
    ! if a gas phase does not exist on either side of the connection, the gas
    ! phase properties from the opposite side are used.
    if (analytical_derivatives) then
      ddelta_pressure_dpup = 1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_up * &
                             wippflo_auxvar_up%d%deng_pg * fmw_comp(iphase)
      ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                             ddensity_kg_ave_dden_kg_dn * &
                             wippflo_auxvar_dn%d%deng_pg * fmw_comp(iphase)
      ddelta_pressure_dpaup = dist_gravity * ddensity_kg_ave_dden_kg_up * &
                              wippflo_auxvar_up%d%deng_pa * fmw_comp(iphase)
      ddelta_pressure_dpadn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                              wippflo_auxvar_dn%d%deng_pa * fmw_comp(iphase)
    endif
    up_scale = 0.d0
    dn_scale = 0.d0
    if (fix_upwind_direction) then
      if (update_upwind_direction .or. count_upwind_direction_flip) then
        prev_upwind_direction = upwind_direction(iphase)
        if (delta_pressure >= 0.d0) then
          ! positive means upstream
          new_upwind_direction = iabs(prev_upwind_direction)
        else
          ! negative means downstream
          new_upwind_direction = -iabs(prev_upwind_direction)
        endif 
        if (count_upwind_direction_flip) then
          if (new_upwind_direction /= prev_upwind_direction) then
            if (derivative_call) then
              gas_upwind_flip_count_by_jac = gas_upwind_flip_count_by_jac + 1
            else
              gas_upwind_flip_count_by_res = gas_upwind_flip_count_by_res + 1
            endif
          endif
        endif
        if (update_upwind_direction) then
          upwind_direction(iphase) = new_upwind_direction
        endif
      endif
      upwind = (upwind_direction(iphase) > 0)
    else
      upwind = (delta_pressure >= 0.d0)
    endif
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
      air_mole_flux = tot_mole_flux
      Res(air_comp_id) = Res(air_comp_id) + air_mole_flux

      if (analytical_derivatives) then
      
        Jgup = 0.d0
        Jgdn = 0.d0
        
        ! Upstream Cell
        ! derivative wrt gas pressure
        ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
        !   liquid pressure derivatives.
        ! derivative total mole flux wrt gas pressure
        dtot_mole_flux_dp = &
          ! ave. liquid density
          q * ddensity_ave_dden_up * wippflo_auxvar_up%d%deng_pg + &
          ! liquid mobility
          up_scale * &
          tot_mole_flux / mobility * wippflo_auxvar_up%d%mobilityg_pg + &
          ! pressure gradient
          tot_mole_flux_ddel_pressure * ddelta_pressure_dpup
        ! derivative water wrt gas pressure
        ! derivative air wrt gas pressure
        Jgup(2,1) = dtot_mole_flux_dp
            
        ! derivative wrt gas saturation
        ! derivative total mole flux wrt gas saturation
        dtot_mole_flux_dsatg = &
          ! relative permeability
          up_scale * &
          tot_mole_flux / mobility * wippflo_auxvar_up%d%mobilityg_satg
        ! derivative water wrt gas saturation
        ! derivative air wrt gas saturation
        Jgup(2,2) = dtot_mole_flux_dsatg
          
        ! Downstream Cell
        ! derivative wrt gas pressure
        ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
        !   liquid pressure derivatives.
        ! derivative total mole flux wrt gas pressure
        dtot_mole_flux_dp = &
          ! ave. liquid density
          q * ddensity_ave_dden_dn * wippflo_auxvar_dn%d%deng_pg + &
          ! liquid mobility
          dn_scale * &
          tot_mole_flux / mobility * wippflo_auxvar_dn%d%mobilityg_pg + &
          ! pressure gradient
          tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
        ! derivative water wrt gas pressure
        ! derivative air wrt gas pressure
        Jgdn(2,1) = dtot_mole_flux_dp
            
        ! derivative wrt gas saturation
        ! derivative total mole flux wrt gas saturation
        dtot_mole_flux_dsatg = &
          ! relative permeability
          dn_scale * &
          tot_mole_flux / mobility * wippflo_auxvar_dn%d%mobilityg_satg
        ! derivative water wrt gas saturation
        ! derivative air wrt gas saturation
        Jgdn(2,2) = dtot_mole_flux_dsatg
          
        Jup = Jup + Jgup
        Jdn = Jdn + Jgdn
      endif
    endif               
  endif

end subroutine WIPPFloFlux

! ************************************************************************** !

subroutine WIPPFloBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         wippflo_auxvar_up,global_auxvar_up, &
                         wippflo_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist,upwind_direction, &
                         wippflo_parameter, &
                         option,v_darcy,Res,J, &
                         derivative_call, &
                         fix_upwind_direction, &
                         update_upwind_direction, &
                         count_upwind_direction_flip, &
                         analytical_derivatives, &
                         debug_connection)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Option_module                              
  use Material_Aux_class
  use General_Aux_module, only : GAS_STATE
  
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
  PetscInt :: upwind_direction(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: J(2,2)
  PetscBool :: derivative_call
  PetscBool :: fix_upwind_direction
  PetscBool :: update_upwind_direction
  PetscBool :: count_upwind_direction_flip
  PetscBool :: analytical_derivatives
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
  PetscInt :: prev_upwind_direction
  PetscInt :: new_upwind_direction
  PetscInt :: iabs_upwind_direction1

  ! Darcy flux
  PetscReal :: ddelta_pressure_dpup, ddelta_pressure_dpdn
  PetscReal :: ddelta_pressure_dpadn
  PetscReal :: dv_darcy_ddelta_pressure
  PetscReal :: dv_darcy_dmobility
  
  PetscReal :: ddensity_kg_ave_dden_kg_up, ddensity_kg_ave_dden_kg_dn
  PetscReal :: ddensity_ave_dden_up, ddensity_ave_dden_dn
  PetscReal :: dtot_mole_flux_dp, dtot_mole_flux_dsatg
  PetscReal :: dpl_dsatg
  PetscReal :: ddelta_pressure_pl
  PetscReal :: tot_mole_flux_ddel_pressure, tot_mole_flux_dmobility
  PetscReal :: dn_scale

  PetscReal :: Jl(2,2)
  PetscReal :: Jg(2,2)
  
  PetscInt :: idof
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id

  Res = 0.d0
  J = 0.d0
  v_darcy = 0.d0  

!geh: we do not want to use the dot product with the unit vector, instead
!     use the principle direction stored in the upwind direction array
!  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  iabs_upwind_direction1 = iabs(upwind_direction(1))
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
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
      if (wippflo_auxvar_up%mobility(iphase) + &
          wippflo_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then
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
        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == GAS_STATE) then
          ! the idea here is to accommodate a free surface boundary
          ! face.  this will not work for an interior grid cell as
          ! there should be capillary pressure in force.
          boundary_pressure = wippflo_auxvar_up%pres(option%gas_phase)
        endif
        density_kg_ave = WIPPFloAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                wippflo_auxvar_up%den_kg, &
                                                wippflo_auxvar_dn%den_kg, &
                                                ddensity_kg_ave_dden_kg_up, &
                                                ddensity_kg_ave_dden_kg_dn)
        ddensity_kg_ave_dden_kg_up = 0.d0 ! always
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          wippflo_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (analytical_derivatives) then
          ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                                 ddensity_kg_ave_dden_kg_dn * &
                                 wippflo_auxvar_dn%d%denl_pl * fmw_comp(iphase)
        endif
        if (bc_type == SEEPAGE_BC .or. &
            bc_type == CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              wippflo_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
            if (analytical_derivatives) then
              option%io_buffer = 'CONDUCTANCE_BC and SEEPAGE_BC need to be &
                &Verified in WIPPFloBCFlux().'
              call printErrMsg(option)
              ddelta_pressure_dpdn = 0.d0
            endif
          endif
        endif
        dn_scale = 0.d0

        if (fix_upwind_direction) then
          if (update_upwind_direction .or. count_upwind_direction_flip) then
            prev_upwind_direction = upwind_direction(iphase)
            if (delta_pressure >= 0.d0) then
              ! positive means upstream
              new_upwind_direction = iabs(prev_upwind_direction)
            else
              ! negative means downstream
              new_upwind_direction = -iabs(prev_upwind_direction)
            endif 
            if (count_upwind_direction_flip) then
              if (new_upwind_direction /= prev_upwind_direction) then
                if (derivative_call) then
                  liq_bc_upwind_flip_count_by_jac = &
                    liq_bc_upwind_flip_count_by_jac + 1
                else
                  liq_bc_upwind_flip_count_by_res = &
                    liq_bc_upwind_flip_count_by_res + 1
                endif
              endif
            endif
            if (update_upwind_direction) then
              upwind_direction(iphase) = new_upwind_direction
            endif
          endif
          upwind = (upwind_direction(iphase) > 0)
        else
          upwind = (delta_pressure >= 0.d0)
        endif
        if (upwind) then
          mobility = wippflo_auxvar_up%mobility(iphase)
        else
          dn_scale = 1.d0        
          mobility = wippflo_auxvar_dn%mobility(iphase)
        endif      

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = WIPPFloAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            wippflo_auxvar_up%den, &
                                            wippflo_auxvar_dn%den, &
                                            ddensity_ave_dden_up, &
                                            ddensity_ave_dden_dn)    
        ddensity_ave_dden_up = 0.d0 ! always
        dv_darcy_dmobility = perm_ave_over_dist * delta_pressure
      endif
    case(NEUMANN_BC)
      dv_darcy_ddelta_pressure = 0.d0
      dv_darcy_dmobility = 0.d0
      ddensity_ave_dden_up = 0.d0
      ddensity_ave_dden_dn = 0.d0
      ddelta_pressure_dpdn = 0.d0
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
          ddensity_ave_dden_dn = 1.d0
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in WIPPFloBCFlux phase loop.'
      call printErrMsg(option)
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
    wat_mole_flux = tot_mole_flux
    Res(wat_comp_id) = Res(wat_comp_id) + wat_mole_flux
   
    if (analytical_derivatives) then
      Jl = 0.d0

      ! derivative wrt gas pressure
      ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
      !   liquid pressure derivatives.
      ! derivative total mole flux wrt gas pressure
      dtot_mole_flux_dp = &
        ! ave. liquid density
        q * ddensity_ave_dden_dn *wippflo_auxvar_dn%d%denl_pl + &
        ! liquid mobility
        dn_scale * &
        tot_mole_flux_dmobility * wippflo_auxvar_dn%d%mobilityl_pl + &
        ! pressure gradient
        tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
      ! derivative water wrt gas pressure
      Jl(1,1) = dtot_mole_flux_dp
      ! derivative air wrt gas pressure
            
      ! derivative wrt gas saturation
      ! pl = pg - pc(satg)
      dpl_dsatg = -1.d0 * wippflo_auxvar_dn%d%pc_satg
      ! delta pressure = plup - pldn
      ddelta_pressure_pl = -1.d0
      ! derivative total mole flux wrt gas saturation
      dtot_mole_flux_dsatg = &
        ! liquid viscosity
        ! since liquid viscosity in a two phase state is a function
        ! of total pressure (gas pressure), there is no derivative
        ! wrt gas saturation
        !dn_scale * &
        !tot_mole_flux_dmobility * &
        !wippflo_auxvar_dn%d%mobilityl_pl * dpl_dsatg + &
        ! relative permeability
        dn_scale * &
        tot_mole_flux_dmobility * &
        wippflo_auxvar_dn%d%mobilityl_satg + &
        !pressure gradient
        tot_mole_flux_ddel_pressure * ddelta_pressure_pl * dpl_dsatg
      ! derivative water wrt gas saturation
      Jl(1,2) = dtot_mole_flux_dsatg
      ! derivative air wrt gas saturation
          
      J = J + Jl
    endif
  endif                   

  iphase = GAS_PHASE
  mobility = 0.d0
  bc_type = ibndtype(iphase)
  select case(bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
      if (wippflo_auxvar_up%mobility(iphase) + &
          wippflo_auxvar_dn%mobility(iphase) > eps) then

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then
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
        if (iphase == LIQUID_PHASE .and. &
            global_auxvar_up%istate == GAS_STATE) then
          ! the idea here is to accommodate a free surface boundary
          ! face.  this will not work for an interior grid cell as
          ! there should be capillary pressure in force.
          boundary_pressure = wippflo_auxvar_up%pres(option%gas_phase)
        endif
        density_kg_ave = WIPPFloAverageDensity(iphase, &
                                                global_auxvar_up%istate, &
                                                global_auxvar_dn%istate, &
                                                wippflo_auxvar_up%den_kg, &
                                                wippflo_auxvar_dn%den_kg, &
                                                ddensity_kg_ave_dden_kg_up, &
                                                ddensity_kg_ave_dden_kg_dn)
        ddensity_kg_ave_dden_kg_up = 0.d0 ! always
        gravity_term = density_kg_ave * dist_gravity
        delta_pressure = boundary_pressure - &
                          wippflo_auxvar_dn%pres(iphase) + &
                          gravity_term
        if (analytical_derivatives) then
          ddelta_pressure_dpadn = dist_gravity * ddensity_kg_ave_dden_kg_dn * &
                                  wippflo_auxvar_dn%d%deng_pa * fmw_comp(iphase)
          ddelta_pressure_dpdn = -1.d0 + dist_gravity * &
                                 ddensity_kg_ave_dden_kg_dn * &
                                 wippflo_auxvar_dn%d%deng_pg * fmw_comp(iphase)
        endif
        if (bc_type == SEEPAGE_BC .or. &
            bc_type == CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              wippflo_auxvar_up%pres(iphase) - &
                option%reference_pressure < eps) then
            delta_pressure = 0.d0
            if (analytical_derivatives) then
              ddelta_pressure_dpdn = 0.d0
            endif
          endif
        endif
        dn_scale = 0.d0
        ! don't expect the derivative to match precisely at delta_pressure = 0
        ! due to potential switch in direction for numerically perturbed
        ! residual
        if (fix_upwind_direction) then
          if (update_upwind_direction .or. count_upwind_direction_flip) then
            prev_upwind_direction = upwind_direction(iphase)
            if (delta_pressure >= 0.d0) then
              ! positive means upstream
              new_upwind_direction = iabs(prev_upwind_direction)
            else
              ! negative means downstream
              new_upwind_direction = -iabs(prev_upwind_direction)
            endif 
            if (count_upwind_direction_flip) then
              if (new_upwind_direction /= prev_upwind_direction) then
                if (derivative_call) then
                 gas_bc_upwind_flip_count_by_jac = &
                   gas_bc_upwind_flip_count_by_jac + 1
                else
                 gas_bc_upwind_flip_count_by_res = &
                   gas_bc_upwind_flip_count_by_res + 1
                endif
              endif
            endif
            if (update_upwind_direction) then
              upwind_direction(iphase) = new_upwind_direction
            endif
          endif
          upwind = (upwind_direction(iphase) > 0)
        else
          upwind = (delta_pressure >= 0.d0)
        endif
        if (upwind) then
          mobility = wippflo_auxvar_up%mobility(iphase)
        else
          dn_scale = 1.d0        
          mobility = wippflo_auxvar_dn%mobility(iphase)
        endif      
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        dv_darcy_ddelta_pressure = perm_ave_over_dist * mobility
        v_darcy(iphase) = dv_darcy_ddelta_pressure * delta_pressure
        ! only need average density if velocity > 0.
        density_ave = WIPPFloAverageDensity(iphase, &
                                            global_auxvar_up%istate, &
                                            global_auxvar_dn%istate, &
                                            wippflo_auxvar_up%den, &
                                            wippflo_auxvar_dn%den, &
                                            ddensity_ave_dden_up, &
                                            ddensity_ave_dden_dn)    
        ddensity_ave_dden_up = 0.d0 ! always
        dv_darcy_dmobility = perm_ave_over_dist * delta_pressure
      endif
    case(NEUMANN_BC)
      dv_darcy_ddelta_pressure = 0.d0
      dv_darcy_dmobility = 0.d0
      ddensity_ave_dden_up = 0.d0 ! always
      ddensity_ave_dden_dn = 0.d0
      ddelta_pressure_dpdn = 0.d0
      ddelta_pressure_dpadn = 0.d0
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
          ddensity_ave_dden_dn = 1.d0
        endif 
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in WIPPFloBCFlux phase loop.'
      call printErrMsg(option)
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
    air_mole_flux = tot_mole_flux
    Res(air_comp_id) = Res(air_comp_id) + air_mole_flux
      
    if (analytical_derivatives) then
      Jg = 0.d0

      ! derivative wrt gas pressure
      ! pl = pg - pc and dpl_dpg = 1.  Therefore, we can use all the 
      !   liquid pressure derivatives.
      ! derivative total mole flux wrt gas pressure
      dtot_mole_flux_dp = &
        ! ave. liquid density
        q * ddensity_ave_dden_dn * wippflo_auxvar_dn%d%deng_pg + &
        ! liquid mobility
        dn_scale * &
        tot_mole_flux_dmobility * wippflo_auxvar_dn%d%mobilityg_pg + &
        ! pressure gradient
        tot_mole_flux_ddel_pressure * ddelta_pressure_dpdn
      ! derivative water wrt gas pressure
      ! derivative air wrt gas pressure
      Jg(2,1) = dtot_mole_flux_dp
            
      ! derivative wrt gas saturation
      ! derivative total mole flux wrt gas saturation
      dtot_mole_flux_dsatg = &
        ! relative permeability
        dn_scale * &
        tot_mole_flux_dmobility * wippflo_auxvar_dn%d%mobilityg_satg
      ! derivative water wrt gas saturation
      ! derivative air wrt gas saturation
      Jg(2,2) = dtot_mole_flux_dsatg
          
      J = J + Jg
    endif
  endif                   

end subroutine WIPPFloBCFlux

! ************************************************************************** !

subroutine WIPPFloSrcSink(option,qsrc,flow_src_sink_type, &
                          wippflo_auxvar,global_auxvar,ss_flow_vol_flux, &
                          scale,Res,J,analytical_derivatives,debug_cell)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module
  
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(wippflo_auxvar_type) :: wippflo_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale
  PetscReal :: Res(option%nflowdof)
  PetscReal :: J(option%nflowdof,option%nflowdof)  
  PetscBool :: analytical_derivatives  
  PetscBool :: debug_cell
      
  PetscReal :: qsrc_mol
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: wat_comp_id, air_comp_id, energy_id
  PetscReal :: Jl(option%nflowdof,option%nflowdof)  
  PetscReal :: Jg(option%nflowdof,option%nflowdof)  
  PetscReal :: dden_bool
  PetscErrorCode :: ierr

  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  
  Res = 0.d0
  J = 0.d0
  
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
  if (analytical_derivatives) then
    Jl = 0.d0
    ! derivative wrt gas pressure
    Jl(1,1) = dden_bool * qsrc(wat_comp_id) * wippflo_auxvar%d%denl_pl
    ! derivative wrt gas saturation
    J = J + Jl
  endif

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
  if (analytical_derivatives) then
    Jg = 0.d0
    ! derivative wrt gas pressure
    Jg(2,1) = dden_bool * qsrc(air_comp_id) * wippflo_auxvar%d%deng_pg
    ! derivative wrt gas saturation
    J = J + Jg
  endif

  if (dabs(qsrc(TWO_INTEGER)) < 1.d-40 .and. &
      qsrc(ONE_INTEGER) < 0.d0) then ! extraction only
    Res(TWO_INTEGER) = qsrc_mol
    ss_flow_vol_flux(air_comp_id) = qsrc_mol/wippflo_auxvar%den(TWO_INTEGER)
    if (analytical_derivatives) then
      Jg = 0.d0
      ! derivative wrt gas pressure
      ! derivative wrt gas saturation
      J = J + Jg
    endif
  endif

end subroutine WIPPFloSrcSink

! ************************************************************************** !

subroutine WIPPFloAccumDerivative(wippflo_auxvar,global_auxvar,material_auxvar, &
                                  soil_heat_capacity,option,J)
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
  PetscReal :: jac(option%nflowdof,option%nflowdof)
  PetscReal :: jac_pert(option%nflowdof,option%nflowdof)
  PetscInt :: idof, irow

  call WIPPFloAccumulation(wippflo_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar,soil_heat_capacity,option, &
                           res,jac,wippflo_analytical_derivatives, &
                           PETSC_FALSE)
                           
  if (wippflo_analytical_derivatives) then
    J = jac
  else
    do idof = 1, option%nflowdof
      call WIPPFloAccumulation(wippflo_auxvar(idof), &
                               global_auxvar, &
                               material_auxvar,soil_heat_capacity, &
                               option,res_pert,jac_pert,PETSC_FALSE,PETSC_FALSE)
      do irow = 1, option%nflowdof
        J(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvar(idof)%pert
      enddo !irow
    enddo ! idof
  endif

end subroutine WIPPFloAccumDerivative

! ************************************************************************** !

subroutine WIPPFloFluxDerivative(wippflo_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 wippflo_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 area, dist, &
                                 upwind_direction, &
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
  
  implicit none
  
  type(wippflo_auxvar_type) :: wippflo_auxvar_up(0:), wippflo_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  PetscInt :: upwind_direction(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_up(option%nflowdof,option%nflowdof)
  PetscReal :: Janal_dn(option%nflowdof,option%nflowdof)
  PetscReal :: Jdummy(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
  option%iflag = -2
  call WIPPFloFlux(wippflo_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   wippflo_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   area,dist,upwind_direction, &
                   wippflo_parameter, &
                   option,v_darcy,res,Janal_up,Janal_dn,&
                   PETSC_TRUE, & ! derivative call 
                   wippflo_fix_upwind_direction, &
                   PETSC_FALSE, & ! update the upwind direction
                   PETSC_FALSE, & ! count upwind direction flip
                   wippflo_analytical_derivatives,PETSC_FALSE)
 
  if (wippflo_analytical_derivatives) then
    Jup = Janal_up
    Jdn = Janal_dn
  else
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call WIPPFloFlux(wippflo_auxvar_up(idof),global_auxvar_up, &
                       material_auxvar_up, &
                       wippflo_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,upwind_direction, &
                       wippflo_parameter, &
                       option,v_darcy,res_pert,Jdummy,Jdummy, &
                       PETSC_TRUE, & ! derivative call
                       wippflo_fix_upwind_direction, &
                       PETSC_FALSE, & ! update the upwind direction
                       wippflo_count_upwind_dir_flip, &
                       PETSC_FALSE,PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jup(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvar_up(idof)%pert
      enddo !irow
    enddo ! idof

    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call WIPPFloFlux(wippflo_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                       material_auxvar_up, &
                       wippflo_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,upwind_direction, &
                       wippflo_parameter, &
                       option,v_darcy,res_pert,Jdummy,Jdummy, &
                       PETSC_TRUE, & ! derivative call
                       wippflo_fix_upwind_direction, &
                       PETSC_FALSE, & ! update the upwind direction
                       wippflo_count_upwind_dir_flip, &
                       PETSC_FALSE,PETSC_FALSE)
      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvar_dn(idof)%pert
      enddo !irow
    enddo ! idof
  endif

end subroutine WIPPFloFluxDerivative

! ************************************************************************** !

subroutine WIPPFloBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   wippflo_auxvar_up, &
                                   global_auxvar_up, &
                                   wippflo_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   area,dist,upwind_direction, &
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
  PetscInt :: upwind_direction(option%nphase)
  type(wippflo_parameter_type) :: wippflo_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow
  PetscReal :: Jdum(option%nflowdof,option%nflowdof)

  Jdn = 0.d0

  option%iflag = -2
  call WIPPFloBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     wippflo_auxvar_up,global_auxvar_up, &
                     wippflo_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist,upwind_direction, &
                     wippflo_parameter, &
                     option,v_darcy,res,Jdum, &
                     PETSC_TRUE, & ! derivative call
                     wippflo_fix_upwind_direction, &
                     PETSC_FALSE, & ! update the upwind direction
                     PETSC_FALSE, & ! count upwind direction flip
                     wippflo_analytical_derivatives,PETSC_FALSE)

  if (wippflo_analytical_derivatives) then
    Jdn = Jdum
  else
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call WIPPFloBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         wippflo_auxvar_up,global_auxvar_up, &
                         wippflo_auxvar_dn(idof),global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist,upwind_direction, &
                         wippflo_parameter, &
                         option,v_darcy,res_pert,Jdum, &
                         PETSC_TRUE, & ! derivative call
                         wippflo_fix_upwind_direction, &
                         PETSC_FALSE, & ! update the upwind direction
                         wippflo_count_upwind_dir_flip, &
                         PETSC_FALSE,PETSC_FALSE)   
      do irow = 1, option%nflowdof
        Jdn(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvar_dn(idof)%pert
      enddo !irow
    enddo ! idof
  endif

end subroutine WIPPFloBCFluxDerivative

! ************************************************************************** !

subroutine WIPPFloSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                    wippflo_auxvars,global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(wippflo_auxvar_type) :: wippflo_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow
  PetscReal :: Jdum(option%nflowdof,option%nflowdof)  

  option%iflag = -3
  call WIPPFloSrcSink(option,qsrc,flow_src_sink_type, &
                      wippflo_auxvars(ZERO_INTEGER),global_auxvar,dummy_real, &
                      scale,res,Jdum,wippflo_analytical_derivatives, &
                      PETSC_FALSE)
                      
  if (wippflo_analytical_derivatives) then
    Jac = Jdum
  else                      
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call WIPPFloSrcSink(option,qsrc,flow_src_sink_type, &
                          wippflo_auxvars(idof),global_auxvar,dummy_real, &
                          scale,res_pert,Jdum,PETSC_FALSE,PETSC_FALSE)            
      do irow = 1, option%nflowdof
        Jac(irow,idof) = (res_pert(irow)-res(irow))/wippflo_auxvars(idof)%pert
      enddo !irow
    enddo ! idof
  endif
  
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

! ************************************************************************** !

subroutine WIPPFloAuxVarDiff(idof,wippflo_auxvar,global_auxvar, &
                             material_auxvar, &
                             wippflo_auxvar_pert,global_auxvar_pert, &
                             material_auxvar_pert, &
                             pert,string,compare_analytical_derivative, &
                             option)

  use Option_module
  use WIPP_Flow_Aux_module
  use Global_Aux_module
  use Material_Aux_class  

  implicit none
  
  type(option_type) :: option
  PetscInt :: idof
  type(wippflo_auxvar_type) :: wippflo_auxvar, wippflo_auxvar_pert
  type(global_auxvar_type) :: global_auxvar, global_auxvar_pert
  class(material_auxvar_type) :: material_auxvar, material_auxvar_pert
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: compare_analytical_derivative
  PetscReal :: pert

  
  PetscInt :: cpid, spid
  PetscInt :: gid, lid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_saturation, gas_saturation
  PetscReal :: liquid_mass_pert, gas_mass_pert
  PetscReal :: liquid_density_pert, gas_density_pert
  PetscReal :: liquid_saturation_pert, gas_saturation_pert
  
  PetscReal :: dpl 
  PetscReal :: dpg 
  PetscReal :: dpc 
  PetscReal :: dps 
  PetscReal :: dsatl
  PetscReal :: dsatg
  PetscReal :: ddenl  
  PetscReal :: ddeng  
  PetscReal :: ddenlkg
  PetscReal :: ddengkg
  PetscReal :: dpsat  
  PetscReal :: dmobilityl  
  PetscReal :: dmobilityg  
  PetscReal :: dmug
  
  PetscReal, parameter :: uninitialized_value = -999.d0
  
  dpl = uninitialized_value
  dpg = uninitialized_value
  dpc = uninitialized_value
  dps = uninitialized_value
  dsatl = uninitialized_value
  dsatg = uninitialized_value
  ddenl = uninitialized_value
  ddeng = uninitialized_value
  ddenlkg = uninitialized_value
  ddengkg = uninitialized_value
  dpsat = uninitialized_value
  dmobilityl = uninitialized_value
  dmobilityg = uninitialized_value
  dmug = uninitialized_value

  lid = option%liquid_phase
  gid = option%gas_phase
  cpid = option%capillary_pressure_id
  spid = option%saturation_pressure_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0
    
  if (compare_analytical_derivative) then
    select case(idof)
      case(1) ! gas pressure pg
        dpl = 1.d0 ! pl = pg - pc
        dpg = 1.d0 ! pg = pg
        dps = 0.d0
        dsatl = 0.d0
        dsatg = 0.d0
        ddenl = wippflo_auxvar%d%denl_pl*dpl
        ddeng = wippflo_auxvar%d%deng_pg
        ddenlkg = ddenl*fmw_comp(1)
        ddengkg = wippflo_auxvar%d%dengkg_pg
            
        dmug = wippflo_auxvar%d%mug_pg
        dmobilityl = wippflo_auxvar%d%mobilityl_pl
        dmobilityg = wippflo_auxvar%d%mobilityg_pg
      case(2) ! gas saturation
        dpl = -1.d0*wippflo_auxvar%d%pc_satg ! pl = pg - pc
        dpg = 0.d0
        dpc = wippflo_auxvar%d%pc_satg
        dps = 0.d0
        dsatl = -1.d0
        dsatg = 1.d0
        ddenl = 0.d0
        dmobilityl = wippflo_auxvar%d%mobilityl_satg
        dmobilityg = wippflo_auxvar%d%mobilityg_satg
    end select
  endif

  print *, '--------------------------------------------------------'
  print *, 'Derivative with respect to ' // trim(string)
  liquid_density = wippflo_auxvar%den(lid)
  gas_density = wippflo_auxvar%den(gid)
  liquid_saturation = wippflo_auxvar%sat(lid)
  gas_saturation = wippflo_auxvar%sat(gid)
  liquid_mass = (liquid_density*liquid_saturation)* &
                 wippflo_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (gas_density*gas_saturation)* & 
              wippflo_auxvar%effective_porosity*material_auxvar%volume
  liquid_density_pert = wippflo_auxvar_pert%den(lid)
  gas_density_pert = wippflo_auxvar_pert%den(gid)
  liquid_saturation_pert = wippflo_auxvar_pert%sat(lid)
  gas_saturation_pert = wippflo_auxvar_pert%sat(gid)
  liquid_mass_pert = (liquid_density_pert*liquid_saturation_pert)* &
              wippflo_auxvar_pert%effective_porosity*material_auxvar_pert%volume
  gas_mass_pert = (gas_density_pert*gas_saturation_pert)* & 
              wippflo_auxvar_pert%effective_porosity*material_auxvar_pert%volume 
  call WIPPFloAuxVarPrintResult('tot liq comp mass [kmol]', &
                                (liquid_mass_pert-liquid_mass)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('tot gas comp mass [kmol]', &
                                (gas_mass_pert-gas_mass)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('         liquid pressure', &
                                (wippflo_auxvar_pert%pres(lid)-wippflo_auxvar%pres(lid))/pert, &
                                dpl,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('            gas pressure', &
                                (wippflo_auxvar_pert%pres(gid)-wippflo_auxvar%pres(gid))/pert, &
                                dpg,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('      capillary pressure', &
                                (wippflo_auxvar_pert%pres(cpid)-wippflo_auxvar%pres(cpid))/pert, &
                                dpc,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('     saturation pressure', &
                                (wippflo_auxvar_pert%pres(spid)-wippflo_auxvar%pres(spid))/pert, &
                                dps,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('       liquid saturation', &
                                (wippflo_auxvar_pert%sat(lid)-wippflo_auxvar%sat(lid))/pert, &
                                dsatl,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('          gas saturation', &
                                (wippflo_auxvar_pert%sat(gid)-wippflo_auxvar%sat(gid))/pert, &
                                dsatg,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('   liquid density [kmol]', &
                                (wippflo_auxvar_pert%den(lid)-wippflo_auxvar%den(lid))/pert, &
                                ddenl,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('      gas density [kmol]', &
                                (wippflo_auxvar_pert%den(gid)-wippflo_auxvar%den(gid))/pert, &
                                ddeng,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('     liquid density [kg]', &
                                (wippflo_auxvar_pert%den_kg(lid)-wippflo_auxvar%den_kg(lid))/pert, &
                                ddenlkg,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('        gas density [kg]', &
                                (wippflo_auxvar_pert%den_kg(gid)-wippflo_auxvar%den_kg(gid))/pert, &
                                ddengkg,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('         temperature [C]', &
                                (wippflo_auxvar_pert%temp-wippflo_auxvar%temp)/pert, &
                                uninitialized_value,uninitialized_value,option)
  !------------------------------
  call WIPPFloAuxVarPrintResult('         liquid mobility', &
                                (wippflo_auxvar_pert%mobility(lid)-wippflo_auxvar%mobility(lid))/pert, &
                                dmobilityl,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('            gas mobility', &
                                (wippflo_auxvar_pert%mobility(gid)-wippflo_auxvar%mobility(gid))/pert, &
                                dmobilityg,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('           gas viscosity', &
                                (wippflo_auxvar_pert%d%mug-wippflo_auxvar%d%mug)/pert, &
                                dmug,uninitialized_value,option)
  call WIPPFloAuxVarPrintResult('      effective porosity', &
                                (wippflo_auxvar_pert%effective_porosity-wippflo_auxvar%effective_porosity)/pert, &
                                uninitialized_value,uninitialized_value,option)
#if 0                                
100 format(a,2(es13.5),es16.8)  
  write(*,100) 'tot liq comp mass [kmol]: ', (liquid_mass_pert-liquid_mass)/pert
  write(*,100) 'tot gas comp mass [kmol]: ', (gas_mass_pert-gas_mass)/pert
  write(*,100) '         liquid pressure: ', (wippflo_auxvar_pert%pres(lid)-wippflo_auxvar%pres(lid))/pert,dpl
  write(*,100) '            gas pressure: ', (wippflo_auxvar_pert%pres(gid)-wippflo_auxvar%pres(gid))/pert,dpg
  write(*,100) '      capillary pressure: ', (wippflo_auxvar_pert%pres(cpid)-wippflo_auxvar%pres(cpid))/pert,dpc
  write(*,100) '     saturation pressure: ', (wippflo_auxvar_pert%pres(spid)-wippflo_auxvar%pres(spid))/pert,dps
  write(*,100) '       liquid saturation: ', (wippflo_auxvar_pert%sat(lid)-wippflo_auxvar%sat(lid))/pert,dsatl
  write(*,100) '          gas saturation: ', (wippflo_auxvar_pert%sat(gid)-wippflo_auxvar%sat(gid))/pert,dsatg
  write(*,100) '   liquid density [kmol]: ', (wippflo_auxvar_pert%den(lid)-wippflo_auxvar%den(lid))/pert,ddenl
  write(*,100) '      gas density [kmol]: ', (wippflo_auxvar_pert%den(gid)-wippflo_auxvar%den(gid))/pert,ddeng
  write(*,100) '     liquid density [kg]: ', (wippflo_auxvar_pert%den_kg(lid)-wippflo_auxvar%den_kg(lid))/pert,ddenl*fmw_comp(1)
  write(*,100) '        gas density [kg]: ', (wippflo_auxvar_pert%den_kg(gid)-wippflo_auxvar%den_kg(gid))/pert,ddengkg
  write(*,100) '         temperature [C]: ', (wippflo_auxvar_pert%temp-wippflo_auxvar%temp)/pert

  write(*,100) '         liquid mobility: ', (wippflo_auxvar_pert%mobility(lid)-wippflo_auxvar%mobility(lid))/pert,dmobilityl
  write(*,100) '            gas mobility: ', (wippflo_auxvar_pert%mobility(gid)-wippflo_auxvar%mobility(gid))/pert,dmobilityg
  write(*,100) '      effective porosity: ', (wippflo_auxvar_pert%effective_porosity-wippflo_auxvar%effective_porosity)/pert
#endif
  write(*,*) '--------------------------------------------------------'  
  
end subroutine WIPPFloAuxVarDiff

! ************************************************************************** !

subroutine WIPPFloAuxVarPrintResult(string,numerical,analytical, &
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
  

end subroutine WIPPFloAuxVarPrintResult

! ************************************************************************** !

subroutine WIPPFloDiffJacobian(string,numerical_jacobian,analytical_jacobian, &
                               residual,residual_pert,perturbation, &
                               perturbation_tolerance,wippflo_auxvar,option)

  use Option_module
  use Utility_module
  
  implicit none
  
  character(len=*) :: string
  PetscReal :: numerical_jacobian(2,2)
  PetscReal :: analytical_jacobian(2,2)
  PetscReal :: residual(2)
  PetscReal :: residual_pert(2,2)
  PetscReal :: perturbation(2)
  PetscReal :: perturbation_tolerance
  type(wippflo_auxvar_type) :: wippflo_auxvar(0:)
  type(option_type) :: option
  
  PetscInt :: irow, icol
  
100 format(2i2,2es13.5,x,a2,es16.8)

  if (len_trim(string) > 1) then
    write(*,'(x,a)') string
  endif
  write(*,'(" Perturbation tolerance: ",es12.4)') perturbation_tolerance
  write(*,'(" r c  numerical    analytical   digits of accuracy")')
  do icol = 1, 2
    do irow = 1, 2
      write(*,100) irow, icol, numerical_jacobian(irow,icol), &
                   analytical_jacobian(irow,icol), &
                   DigitsOfAccuracy(numerical_jacobian(irow,icol), &
                                    analytical_jacobian(irow,icol))
    enddo
  enddo

#if 0
200 format(2es20.12)
300 format(a24,10es20.12)
  do icol = 1, 2
    write(*,'(/," dof = ",i1,"  perturbation = ",es13.5)') icol, perturbation(icol)
!    write(*,300) 'density', wippflo_auxvar(icol)%den(:), wippflo_auxvar(0)%den(:)
!    write(*,300) 'energy', wippflo_auxvar(icol)%U(:), wippflo_auxvar(0)%U(:)
    write(*,'("  residual_pert       residual")')
    do irow = 1, 2
      write(*,200) residual_pert(irow,icol), residual(irow)
    enddo
  enddo
#endif  
  
end subroutine WIPPFloDiffJacobian

end module WIPP_Flow_Common_module
