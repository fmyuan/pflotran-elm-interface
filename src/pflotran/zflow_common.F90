module ZFlow_Common_module

  use ZFlow_Aux_module
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
    subroutine FluxDummy(zflow_auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         zflow_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area, dist, &
                         zflow_parameter, &
                         option,v_darcy,Res,Jup,Jdn, &
                         dResdparamup,dResdparamdn, &
                         calculate_derivatives)
      use ZFlow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(zflow_parameter_type) :: zflow_parameter
      PetscReal :: Res(ZFLOW_MAX_DOF)
      PetscReal :: Jup(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
      PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
      PetscReal :: dResdparamup(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
      PetscReal ::  dResdparamdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
      PetscBool :: calculate_derivatives
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                           zflow_auxvar_up,global_auxvar_up, &
                           zflow_auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           area,dist, &
                           zflow_parameter, &
                           option,v_darcy,Res,Jdn, &
                           dResdparamdn, &
                           calculate_derivatives)
      use ZFlow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(option_type) :: option
      PetscInt :: ibndtype(1)
      PetscInt :: auxvar_mapping(ZFLOW_MAX_INDEX)
      PetscReal :: auxvars(:) ! from aux_real_var array
      type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_dn
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(zflow_parameter_type) :: zflow_parameter
      PetscReal :: v_darcy
      PetscReal :: Res(ZFLOW_MAX_DOF)
      PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
      PetscReal :: dResdparamdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
      PetscBool :: calculate_derivatives
    end subroutine BCFluxDummy
  end interface

  public :: ZFlowAccumulation, &
            ZFlowFluxHarmonicPermOnly, &
            ZFlowBCFluxHarmonicPermOnly, &
            ZFlowAccumDerivative, &
            XXFluxDerivative, &
            XXBCFluxDerivative, &
            ZFlowSrcSinkDerivative

  public :: XXFlux, &
            XXBCFlux

contains

! ************************************************************************** !

subroutine ZFlowAccumulation(zflow_auxvar,global_auxvar,material_auxvar, &
                             option,Res,Jac,dResdparam,calculate_derivatives)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jac(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparam(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: porosity
  PetscReal :: saturation
  PetscReal :: por_sat
  PetscReal :: volume_over_dt
  PetscReal :: tempreal

  Res = 0.d0
  Jac = 0.d0
  dResdparam = 0.d0

    ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use zflow_auxvar%effective porosity here as it enables numerical
  ! derivatives to be employed
  porosity = zflow_auxvar%effective_porosity
  saturation = zflow_auxvar%sat
  por_sat = saturation * porosity

  if (zflow_liq_flow_eq > 0) then
    ! accumulation term units = m^3 liquid/s
    ! Res[m^3 liquid/sec] = sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
    !                       vol[m^3 bulk] / dt[sec]
    tempreal = volume_over_dt * zflow_density_kmol
    Res(zflow_liq_flow_eq) = por_sat * tempreal

    if (calculate_derivatives) then
      Jac(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
        tempreal * (zflow_auxvar%dsat_dp * porosity + &
                    saturation * zflow_auxvar%dpor_dp)
      if (zflow_calc_adjoint) then
        if (zflow_adjoint_parameter == ZFLOW_ADJOINT_POROSITY) then
          dResdparam(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
            tempreal * zflow_auxvar%sat
        endif
      endif
    endif
  endif

  if (zflow_sol_tran_eq > 0) then
    ! accumulation term units = mol/s
    ! Res[mole/sec] = c [mol/L] * 1000 [L/m^3 liquid]
    !                 sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
    !                 vol[m^3 bulk] / dt[sec]
    tempreal = 1000.d0 * volume_over_dt
    Res(zflow_sol_tran_eq) = zflow_auxvar%conc * tempreal * por_sat

    if (calculate_derivatives) then
      Jac(zflow_sol_tran_eq,zflow_sol_tran_eq) = tempreal * por_sat
      if (zflow_liq_flow_eq > 0) then
        Jac(zflow_sol_tran_eq,zflow_liq_flow_eq) = &
          zflow_auxvar%conc * tempreal * &
          (zflow_auxvar%dsat_dp * porosity + &
          saturation * zflow_auxvar%dpor_dp)
      endif
      if (zflow_calc_adjoint) then
        if (zflow_adjoint_parameter == ZFLOW_ADJOINT_POROSITY) then
          dResdparam(zflow_sol_tran_eq,zflow_sol_tran_eq) = &
            zflow_auxvar%conc * tempreal * zflow_auxvar%sat
        endif
      endif
    endif
  endif

end subroutine ZFlowAccumulation

! ************************************************************************** !

subroutine ZFlowFluxHarmonicPermOnly(zflow_auxvar_up,global_auxvar_up, &
                                     material_auxvar_up, &
                                     zflow_auxvar_dn,global_auxvar_dn, &
                                     material_auxvar_dn, &
                                     area, dist, &
                                     zflow_parameter, &
                                     option,v_darcy,Res, &
                                     Jup,Jdn, &
                                     dResdparamup,dResdparamdn, &
                                     calculate_derivatives)
  !
  ! Computes the internal flux terms for the residual based on harmonic
  ! intrinsic permeability
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !
  use Option_module
  use Material_Aux_module
  use Connection_module

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jup(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparamup(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparamdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscReal :: perm_ave_over_dist_visc
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: dkr_dpup, dkr_dpdn
  PetscReal :: tempreal
  PetscReal :: numerator, denominator
  PetscReal :: dperm_ave_dKup, dperm_ave_dKdn
  PetscReal :: delta_conc
  PetscReal :: dq_dpup, dq_dpdn
  PetscReal :: conc_upwind, dconc_upwind_dup, dconc_upwind_ddn
  PetscReal :: D_hyd_up, D_hyd_dn
  PetscReal, parameter :: D_molecular = 0.d0
  PetscReal, parameter :: D_mech_up = 0.d0, D_mech_dn = 0.d0
  PetscReal :: Deff_over_dist
  PetscReal :: dD_mech_up_dpup, dD_mech_dn_dpdn
  PetscReal :: dDeff_over_dist_dpup, dDeff_over_dist_dpdn
  PetscReal :: dD_hyd_up_dpup, dD_hyd_dn_dpdn

  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  dq_dpup = 0.d0
  dq_dpdn = 0.d0
  dResdparamup = 0.d0
  dResdparamdn = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)

  if (zflow_liq_flow_eq > 0) then

    call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
    call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

    numerator = perm_up * perm_dn
    denominator = dist_up*perm_dn + dist_dn*perm_up
    perm_ave_over_dist_visc = numerator / (denominator * zflow_viscosity)

    if (zflow_auxvar_up%kr + zflow_auxvar_dn%kr > eps) then

      gravity_term = zflow_density_kg * dist_gravity
      delta_pressure = zflow_auxvar_up%pres - &
                       zflow_auxvar_dn%pres + &
                       gravity_term
      dkr_dpup = 0.d0
      dkr_dpdn = 0.d0
      if (zflow_tensorial_rel_perm) then
        if (delta_pressure >= 0.d0) then
          call ZFlowAuxTensorialRelPerm(zflow_auxvar_up, &
                  zflow_parameter% &
                    tensorial_rel_perm_exponent(:,material_auxvar_up%id), &
                  dist,kr,dkr_dpup,option)
        else
          call ZFlowAuxTensorialRelPerm(zflow_auxvar_dn, &
                  zflow_parameter% &
                    tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                  dist,kr,dkr_dpdn,option)
        endif
      else
        if (delta_pressure >= 0.d0) then
          kr = zflow_auxvar_up%kr
          dkr_dpup = zflow_auxvar_up%dkr_dp
        else
          kr = zflow_auxvar_dn%kr
          dkr_dpdn = zflow_auxvar_dn%dkr_dp
        endif
      endif

      if (kr > floweps) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy = perm_ave_over_dist_visc * kr * delta_pressure
        ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy * area
        ! Res[m^3 liquid/sec]
        Res(zflow_liq_flow_eq) = Res(zflow_liq_flow_eq) + &
                                                 q * zflow_density_kmol

        if (calculate_derivatives) then
          tempreal = perm_ave_over_dist_visc * area
          dq_dpup = tempreal * (delta_pressure * dkr_dpup + kr)
          dq_dpdn = tempreal * (delta_pressure * dkr_dpdn - kr)
          Jup(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
                                              dq_dpup * zflow_density_kmol
          Jdn(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
                                              dq_dpdn * zflow_density_kmol
          if (zflow_calc_adjoint) then
            if (zflow_adjoint_parameter == ZFLOW_ADJOINT_PERMEABILITY) then
              tempreal = denominator * denominator * zflow_viscosity
              dperm_ave_dKup = perm_dn * perm_dn * dist_up / tempreal
              dperm_ave_dKdn = perm_up * perm_up * dist_dn / tempreal
              tempreal = kr * delta_pressure * area * zflow_density_kmol
              dResdparamup(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
                tempreal * dperm_ave_dKup
              dResdparamdn(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
                tempreal * dperm_ave_dKdn
            endif
          endif
        endif
      endif
    endif
  endif

  if (zflow_sol_tran_eq > 0) then
    if (q >= 0) then
      conc_upwind = zflow_auxvar_up%conc
      dconc_upwind_dup = 1.d0
      dconc_upwind_ddn = 0.d0
    else
      conc_upwind = zflow_auxvar_dn%conc
      dconc_upwind_dup = 0.d0
      dconc_upwind_ddn = 1.d0
    endif
    delta_conc = zflow_auxvar_up%conc - zflow_auxvar_dn%conc
    D_hyd_up = D_mech_up + &
               zflow_auxvar_up%effective_porosity * zflow_auxvar_up%sat * &
               material_auxvar_up%tortuosity * D_molecular
    D_hyd_dn = D_mech_dn + &
               zflow_auxvar_dn%effective_porosity * zflow_auxvar_dn%sat * &
               material_auxvar_dn%tortuosity * D_molecular
    numerator = D_hyd_up * D_hyd_dn
    denominator = dist_up*D_hyd_dn + D_hyd_up*D_hyd_up
    if (D_hyd_dn > 0.d0) then
      option%io_buffer = 'Update denominator'
      call PrintErrMsg(option)
    else
      denominator = 1.d0
    endif
    Deff_over_dist = numerator / denominator
    ! Res[mol/sec]
    Res(zflow_sol_tran_eq) = Res(zflow_sol_tran_eq) + &
      (q * 1000.d0 * conc_upwind - & ! advection
       area * Deff_over_dist * delta_conc) ! hydrodynamic dispersion
    if (calculate_derivatives) then
      Jup(zflow_sol_tran_eq,zflow_sol_tran_eq) = &
        (q * 1000.d0 * dconc_upwind_dup - &
         area * Deff_over_dist * 1.d0)
      Jdn(zflow_sol_tran_eq,zflow_sol_tran_eq) = &
        (q * 1000.d0 * dconc_upwind_ddn - &
         area * Deff_over_dist * (-1.d0))
      if (zflow_liq_flow_eq > 0) then
        if (D_mech_up > 0.d0) then
          ! implement derivatives here and in boundary flux
          option%io_buffer = 'Implement derivatives for D_mech_up'
          call PrintErrMsg(option)
        endif
        dD_mech_up_dpup = 0.d0
        dD_hyd_up_dpup = dD_mech_up_dpup + &
                         (zflow_auxvar_up%dsat_dp * &
                          zflow_auxvar_up%effective_porosity + &
                          zflow_auxvar_up%sat * zflow_auxvar_up%dpor_dp) * &
                         material_auxvar_up%tortuosity * D_molecular
        dD_mech_dn_dpdn = 0.d0
        dD_hyd_dn_dpdn = dD_mech_dn_dpdn + &
                         (zflow_auxvar_dn%dsat_dp * &
                          zflow_auxvar_dn%effective_porosity + &
                          zflow_auxvar_dn%sat * zflow_auxvar_dn%dpor_dp) * &
                         material_auxvar_dn%tortuosity * D_molecular
        tempreal = denominator * denominator
        dDeff_over_dist_dpup = dD_hyd_up_dpup * D_hyd_dn * D_hyd_dn * &
                               dist_up / tempreal
        dDeff_over_dist_dpdn = dD_hyd_dn_dpdn * D_hyd_up * D_hyd_up * &
                               dist_dn / tempreal
        Jup(zflow_sol_tran_eq,zflow_liq_flow_eq) = &
          (dq_dpup * 1000.d0 * conc_upwind - &
           area * dDeff_over_dist_dpup * 1.d0)
        Jdn(zflow_sol_tran_eq,zflow_liq_flow_eq) = &
          (dq_dpdn * 1000.d0 * conc_upwind - &
           area * dDeff_over_dist_dpdn * (-1.d0))
      endif
    endif
  endif

end subroutine ZFlowFluxHarmonicPermOnly

! ************************************************************************** !

subroutine ZFlowBCFluxHarmonicPermOnly(ibndtype,auxvar_mapping,auxvars, &
                                       zflow_auxvar_up,global_auxvar_up, &
                                       zflow_auxvar_dn,global_auxvar_dn, &
                                       material_auxvar_dn, &
                                       area,dist, &
                                       zflow_parameter, &
                                       option,v_darcy,Res,Jdn, &
                                       dResdparamdn, &
                                       calculate_derivatives)
  !
  ! Computes the boundary flux terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !
  use Option_module
  use Material_Aux_module
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1)
  PetscInt :: auxvar_mapping(ZFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparamdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscInt :: bc_type
  PetscInt :: idof
  PetscReal :: perm_ave_over_dist_visc
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: dkr_dp
  PetscReal :: perm_dn
  PetscReal :: boundary_pressure
  PetscReal :: ddelta_pressure_dpdn
  PetscBool :: derivative_toggle
  PetscReal :: dperm_dK
  PetscReal :: tempreal
  PetscReal :: delta_conc
  PetscReal :: dq_dpdn
  PetscReal :: conc_upwind, dconc_upwind_ddn
  PetscReal :: D_hyd_dn
  PetscReal, parameter :: D_molecular = 0.d0
  PetscReal, parameter :: D_mech_dn = 0.d0
  PetscReal :: Deff_over_dist
  PetscReal :: dD_mech_dn_dpdn
  PetscReal :: dDeff_over_dist_dpdn
  PetscReal :: dD_hyd_dn_dpdn


  Res = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  dq_dpdn = 0.d0
  dperm_dK = 0.d0
  dResdparamdn = 0.d0

  if (zflow_liq_flow_eq > 0) then
    call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

    kr = 0.d0
    bc_type = ibndtype(auxvar_mapping(ZFLOW_COND_WATER_INDEX))
    derivative_toggle = PETSC_TRUE
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
          HYDROSTATIC_CONDUCTANCE_BC)
        if (zflow_auxvar_up%kr + zflow_auxvar_dn%kr > eps) then
          ! dist(0) = scalar - magnitude of distance
          ! gravity = vector(3)
          ! dist(1:3) = vector(3) - unit vector
          dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

          if (bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
            idof = auxvar_mapping(ZFLOW_COND_CONDUCTANCE_INDEX)
            perm_ave_over_dist_visc = auxvars(idof) / zflow_viscosity
          else
            perm_ave_over_dist_visc = perm_dn / (dist(0) * zflow_viscosity)
            dperm_dK = 1.d0 / (dist(0) * zflow_viscosity)
          endif

          boundary_pressure = zflow_auxvar_up%pres
          gravity_term = zflow_density_kg * dist_gravity
          delta_pressure = boundary_pressure - &
                          zflow_auxvar_dn%pres + &
                          gravity_term
          ddelta_pressure_dpdn = -1.d0
          if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
              bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                zflow_auxvar_up%pres - &
                  option%flow%reference_pressure < eps) then
              delta_pressure = 0.d0
              ddelta_pressure_dpdn = 0.d0
            endif
          endif
          if (zflow_tensorial_rel_perm) then
            if (delta_pressure >= 0.d0) then
              call ZFlowAuxTensorialRelPerm(zflow_auxvar_up, &
                      zflow_parameter% &
                        tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                      dist,kr,dkr_dp,option)
              ! override as the upwind pressure is fixed
              dkr_dp = 0.d0
            else
              call ZFlowAuxTensorialRelPerm(zflow_auxvar_dn, &
                      zflow_parameter% &
                        tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                      dist,kr,dkr_dp,option)
            endif
          else
            if (delta_pressure >= 0.d0) then
              kr = zflow_auxvar_up%kr
              dkr_dp = 0.d0
            else
              kr = zflow_auxvar_dn%kr
              dkr_dp = zflow_auxvar_dn%dkr_dp
            endif
          endif

          ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
          !                    dP[Pa]]
          v_darcy = perm_ave_over_dist_visc * kr * delta_pressure
        endif

      case(NEUMANN_BC)
        idof = auxvar_mapping(ZFLOW_COND_WATER_INDEX)
        if (dabs(auxvars(idof)) > floweps) then
          v_darcy = auxvars(idof)
        endif
        derivative_toggle = PETSC_FALSE
      case default
        option%io_buffer = &
          'Boundary condition type (' // trim(StringWrite(bc_type)) // &
          ') not recognized in ZFlowBCFlux phase loop.'
        call PrintErrMsg(option)
    end select
    if (dabs(v_darcy) > 0.d0 .or. kr > 0.d0) then
      ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy * area
      ! Res[m^3 liquid/sec]
      Res(zflow_liq_flow_eq) = Res(zflow_liq_flow_eq) + q * zflow_density_kmol
      if (calculate_derivatives .and. derivative_toggle) then
        ! derivative toggle takes care of the NEUMANN side
        dq_dpdn = area * perm_ave_over_dist_visc * &
                  (dkr_dp * delta_pressure + ddelta_pressure_dpdn * kr)
        Jdn(zflow_liq_flow_eq,zflow_liq_flow_eq) = dq_dpdn * zflow_density_kmol
        if (zflow_calc_adjoint) then
          if (zflow_adjoint_parameter == ZFLOW_ADJOINT_PERMEABILITY) then
            tempreal = area * zflow_density_kmol * kr
            dResdparamdn(zflow_liq_flow_eq,zflow_liq_flow_eq) = &
              tempreal * delta_pressure * dperm_dK
          endif
        endif
      endif
    endif
  endif

  if (zflow_sol_tran_eq > 0) then
    if (q >= 0) then
      conc_upwind = zflow_auxvar_up%conc
      dconc_upwind_ddn = 0.d0
    else
      conc_upwind = zflow_auxvar_dn%conc
      dconc_upwind_ddn = 1.d0
    endif
    delta_conc = zflow_auxvar_up%conc - zflow_auxvar_dn%conc
    D_hyd_dn = D_mech_dn + &
               zflow_auxvar_dn%effective_porosity * zflow_auxvar_dn%sat * &
               material_auxvar_dn%tortuosity * D_molecular
    Deff_over_dist = D_hyd_dn / dist(0)
    ! Res[mol/sec]
    Res(zflow_sol_tran_eq) = Res(zflow_sol_tran_eq) + &
      (q * 1000.d0 * conc_upwind - & ! advection
       area * Deff_over_dist * delta_conc) ! hydrodynamic dispersion
    if (calculate_derivatives) then
      Jdn(zflow_sol_tran_eq,zflow_sol_tran_eq) = &
        (q * 1000.d0 * dconc_upwind_ddn - &
         area * Deff_over_dist * (-1.d0))
      if (zflow_liq_flow_eq > 0) then
        dD_mech_dn_dpdn = 0.d0
        dD_hyd_dn_dpdn = dD_mech_dn_dpdn + &
                         (zflow_auxvar_dn%dsat_dp * &
                          zflow_auxvar_dn%effective_porosity + &
                          zflow_auxvar_dn%sat * zflow_auxvar_dn%dpor_dp) * &
                         material_auxvar_dn%tortuosity * D_molecular
        dDeff_over_dist_dpdn = dD_hyd_dn_dpdn / dist(0)
        Jdn(zflow_sol_tran_eq,zflow_liq_flow_eq) = &
          (dq_dpdn * 1000.d0 * conc_upwind - &
           area * dDeff_over_dist_dpdn * (-1.d0))
      endif
    endif
  endif

end subroutine ZFlowBCFluxHarmonicPermOnly

! ************************************************************************** !

subroutine ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                        zflow_auxvar,global_auxvar,material_auxvar, &
                        ss_flow_vol_flux,scale,Res,Jdn, &
                        calculate_derivatives)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !
  use Option_module
  use Material_Aux_module
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: scale
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: qsrc_m3, qsrc_L

  Res = 0.d0
  Jdn = 0.d0

  if (zflow_liq_flow_eq > 0) then
    ! liquid phase
    qsrc_m3 = 0.d0
    select case(flow_src_sink_type)
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec
        qsrc_m3 = qsrc(zflow_liq_flow_eq)
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_m3 = qsrc(zflow_liq_flow_eq)*scale
      case default
        option%io_buffer = 'src_sink_type not supported in ZFlowSrcSink'
        call PrintErrMsg(option)
    end select
    ss_flow_vol_flux = qsrc_m3
    ! Res[m^3 liquid/sec]
    Res(zflow_liq_flow_eq) = qsrc_m3 * zflow_density_kmol

    if (calculate_derivatives) then
    endif
  endif

  if (zflow_sol_tran_eq > 0) then
    qsrc_L = qsrc_m3 * 1000.d0
    Res(zflow_sol_tran_eq) = qsrc_L * qsrc(zflow_sol_tran_eq)
    if (calculate_derivatives) then
      Jdn(zflow_sol_tran_eq,zflow_sol_tran_eq) = qsrc_L
      if (zflow_liq_flow_eq > 0) then
      endif
    endif
  endif

end subroutine ZFlowSrcSink

! ************************************************************************** !

subroutine ZFlowAccumDerivative(zflow_auxvar,global_auxvar, &
                                material_auxvar, &
                                option,Res,Jac,dResdparam)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jac(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparam(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)

  PetscReal :: res_pert(ZFLOW_MAX_DOF)
  PetscReal :: Jdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dJdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  call ZFlowAccumulation(zflow_auxvar(ZERO_INTEGER), &
                         global_auxvar, &
                         material_auxvar, &
                         option,Res,Jac,dResdparam, &
                         .not.zflow_numerical_derivatives)

  if (zflow_numerical_derivatives) then
    do idof = 1, option%nflowdof
      call ZFlowAccumulation(zflow_auxvar(idof), &
                            global_auxvar, &
                            material_auxvar, &
                            option,res_pert,Jdum,dJdum, &
                            PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq))/zflow_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine ZFlowAccumDerivative

! ************************************************************************** !

subroutine XXFluxDerivative(zflow_auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            zflow_auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            area, dist,zflow_parameter,option,v_darcy, &
                            Res,Jup,Jdn,dResdparamup,dResdparamdn)
  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar_up(0:), zflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jup(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparamup(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparamdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)

  PetscReal :: res_pert(ZFLOW_MAX_DOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dJdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  Jup = 0.d0
  Jdn = 0.d0

  call XXFlux(zflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist,zflow_parameter,option,v_darcy, &
              Res,Jup,Jdn,dResdparamup,dResdparamdn, &
              .not.zflow_numerical_derivatives)

  if (zflow_numerical_derivatives) then
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(zflow_auxvar_up(idof),global_auxvar_up, &
                  material_auxvar_up, &
                  zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,zflow_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jup(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        zflow_auxvar_up(idof)%pert
      enddo
    enddo
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(zflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  zflow_auxvar_dn(idof),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,zflow_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        zflow_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine XXFluxDerivative

! ************************************************************************** !

subroutine XXBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                              zflow_auxvar_up, &
                              global_auxvar_up, &
                              zflow_auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              area,dist,zflow_parameter,option,v_darcy, &
                              Res,Jdn,dResdparamdn)
  !
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:1)
  PetscInt :: auxvar_mapping(ZFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dResdparamdn(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)

  PetscReal :: res_pert(ZFLOW_MAX_DOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscReal :: dJdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  Jdn = 0.d0

  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                zflow_auxvar_up,global_auxvar_up, &
                zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                material_auxvar_dn, &
                area,dist,zflow_parameter,option,v_darcy, &
                Res,Jdn,dResdparamdn,.not.zflow_numerical_derivatives)

  if (zflow_numerical_derivatives) then
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                    zflow_auxvar_up,global_auxvar_up, &
                    zflow_auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    area,dist,zflow_parameter,option,v_darcy_dum, &
                    res_pert,Jdum,dJdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        zflow_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine XXBCFluxDerivative

! ************************************************************************** !

subroutine ZFlowSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                  zflow_auxvars,global_auxvar, &
                                  material_auxvar, &
                                  ss_flow_vol_flux,scale, &
                                  Res,Jac)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 01/10/22
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type
  type(zflow_auxvar_type) :: zflow_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: scale
  PetscReal :: Res(ZFLOW_MAX_DOF)
  PetscReal :: Jac(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)

  PetscReal :: res_pert(ZFLOW_MAX_DOF)
  PetscReal :: dummy_real
  PetscReal :: Jdum(ZFLOW_MAX_DOF,ZFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  Jac = 0.d0
  ! unperturbed zflow_auxvars value
  call ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                    zflow_auxvars(ZERO_INTEGER),global_auxvar, &
                    material_auxvar, &
                    ss_flow_vol_flux,scale, &
                    Res,Jac,.not.zflow_numerical_derivatives)

  if (zflow_numerical_derivatives) then
    ! perturbed zflow_auxvars values
    do idof = 1, option%nflowdof
      call ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                        zflow_auxvars(idof),global_auxvar, &
                        material_auxvar, &
                        dummy_real,scale, &
                        res_pert,Jdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
      Jac(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                      zflow_auxvars(idof)%pert
      enddo
    enddo
  endif

end subroutine ZFlowSrcSinkDerivative

end module ZFlow_Common_module
