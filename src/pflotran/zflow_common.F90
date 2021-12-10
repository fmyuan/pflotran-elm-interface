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
                         dJupdKup,dJupdKdn,dJdndKup,JdndKdn, &
                         drhsdKup,drhsdKdn, &
                         dResdKup,dResdKdn, &
                         calculate_derivatives)
      use ZFlow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      implicit none
      type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy(1)
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(zflow_parameter_type) :: zflow_parameter
      PetscReal :: Res(1)
      PetscReal :: Jup(1,1), Jdn(1,1)
      PetscReal :: dJupdKup(1), dJupdKdn(1), dJdndKup(1), JdndKdn(1)
      PetscReal :: drhsdKup(1), drhsdKdn(1)
      PetscReal :: dResdKup(1), dResdKdn(1)
      PetscBool :: calculate_derivatives
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                           zflow_auxvar_up,global_auxvar_up, &
                           zflow_auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           area,dist, &
                           zflow_parameter, &
                           option,v_darcy,Res,Jdn, &
                           dJdndKdn,drhsdKdn,dResdKdn, &
                           calculate_derivatives)
      use ZFlow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      implicit none
      type(option_type) :: option
      PetscInt :: ibndtype(1)
      PetscInt :: auxvar_mapping(ZFLOW_MAX_INDEX)
      PetscReal :: auxvars(:) ! from aux_real_var array
      type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_dn
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(zflow_parameter_type) :: zflow_parameter
      PetscReal :: v_darcy(1)
      PetscReal :: Res(1)
      PetscReal :: Jdn(1,1)
      PetscReal :: dJdndKdn(1)
      PetscReal :: drhsdKdn(1)
      PetscReal :: dResdKdn(1)
      PetscBool :: calculate_derivatives
    end subroutine BCFluxDummy
  end interface

  public :: ZFlowAccumulation, &
            ZFlowFluxHarmonicPermOnly, &
            ZFlowBCFluxHarmonicPermOnly, &
            ZFlowSrcSink, &
            ZFlowAccumDerivative, &
            XXFluxDerivative, &
            XXBCFluxDerivative, &
            ZFlowSrcSinkDerivative

  public :: XXFlux, &
            XXBCFlux

contains

! ************************************************************************** !

subroutine ZFlowAccumulation(zflow_auxvar,global_auxvar,material_auxvar, &
                             option,Res,Jac,dResdpor,calculate_derivatives)
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

  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(1)
  PetscReal :: Jac(1,1)
  PetscReal :: dResdpor(1)
  PetscBool :: calculate_derivatives

  PetscReal :: porosity
  PetscReal :: volume_over_dt

  Res = 0.d0
  Jac = 0.d0
  dResdpor = 0.d0

  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use zflow_auxvar%effective porosity here as it enables numerical
  ! derivatives to be employed
  porosity = zflow_auxvar%effective_porosity

  ! accumulation term units = m^3 liquid/s
  ! Res[m^3 liquid/sec] = sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
  !                       vol[m^3 bulk] / dt[sec]
  Res(1) = zflow_auxvar%sat * porosity * volume_over_dt * zflow_density_kmol

  if (calculate_derivatives) then
    Jac(1,1) = Res(1) / zflow_auxvar%sat * zflow_auxvar%dsat_dp + &
               Res(1) / porosity * zflow_auxvar%dpor_dp
    if (zflow_calc_adjoint) then
      if (zflow_adjoint_parameter == ZFLOW_ADJOINT_POROSITY) then
        dResdpor(1) = Res(1)/porosity
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
                                     dJupdKup,dJupdKdn,dJdndKup,dJdndKdn, &
                                     drhsdKup,drhsdKdn, &
                                     dResdKup,dResdKdn, &
                                     calculate_derivatives)
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

  implicit none

  type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(1)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: Res(1)
  PetscReal :: Jup(1,1), Jdn(1,1)
  PetscReal :: dJupdKup(1), dJupdKdn(1), dJdndKup(1), dJdndKdn(1)
  PetscReal :: drhsdKup(1), drhsdKdn(1)
  PetscReal :: dResdKup(1), dResdKdn(1)
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

  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)


  numerator = perm_up * perm_dn
  denominator = dist_up*perm_dn + dist_dn*perm_up
  perm_ave_over_dist_visc = numerator / (denominator * zflow_viscosity)

  if (zflow_auxvar_up%kr + zflow_auxvar_dn%kr > eps) then

    gravity_term = zflow_density_kg * dist_gravity
    delta_pressure = zflow_auxvar_up%pres - &
                     zflow_auxvar_dn%pres + &
                     gravity_term
    if (delta_pressure >= 0.d0) then
      kr = zflow_auxvar_up%kr
      dkr_dpup = zflow_auxvar_up%dkr_dp
      dkr_dpdn = 0.d0
    else
      kr = zflow_auxvar_dn%kr
      dkr_dpup = 0.d0
      dkr_dpdn = zflow_auxvar_dn%dkr_dp
    endif

    if (kr > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(1) = perm_ave_over_dist_visc * kr * delta_pressure
      ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(1) * area * zflow_density_kmol
      ! Res[m^3 liquid/sec]
      Res = Res + q

      if (calculate_derivatives) then
        tempreal = perm_ave_over_dist_visc * area * zflow_density_kmol
        Jup(1,1) = tempreal * (delta_pressure * dkr_dpup + kr)
        Jdn(1,1) = tempreal * (delta_pressure * dkr_dpdn - kr)
        if (zflow_calc_adjoint) then
          if (zflow_adjoint_parameter == ZFLOW_ADJOINT_PERMEABILITY) then
            tempreal = denominator * denominator * zflow_viscosity
            dperm_ave_dKup = perm_dn * perm_dn * dist_up / tempreal
            dperm_ave_dKdn = perm_up * perm_up * dist_dn / tempreal
            tempreal = area * zflow_density_kmol * kr
            dJupdKup(1) = dperm_ave_dKup * tempreal
            dJupdKdn(1) = dperm_ave_dKdn * tempreal
            tempreal = -1.d0 * area * zflow_density_kmol * kr
            dJdndKup(1) = dperm_ave_dKup * tempreal
            dJdndKdn(1) = dperm_ave_dKdn * tempreal
            tempreal = area * zflow_density_kmol * kr * gravity_term
            drhsdKup(1) = dperm_ave_dKup * tempreal
            drhsdKdn(1) = dperm_ave_dKdn * tempreal
            tempreal = Res(1) / perm_ave_over_dist_visc
            dResdKup(1) = tempreal * dperm_ave_dKup
            dResdKdn(1) = tempreal * dperm_ave_dKdn
          endif
        endif
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
                                       dJdndKdn,drhsdKdn,dResdKdn, &
                                       calculate_derivatives)
  !
  ! Computes the boundary flux terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  !
  use Option_module
  use Material_Aux_class
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1)
  PetscInt :: auxvar_mapping(ZFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: v_darcy(1)
  PetscReal :: Res(1)
  PetscReal :: Jdn(1,1)
  PetscReal :: dJdndKdn(1)
  PetscReal :: drhsdKdn(1)
  PetscReal :: dResdKdn(1)
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

  Res = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0

  dperm_dK = 0.d0

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  kr = 0.d0
  bc_type = ibndtype(1)
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
          idof = auxvar_mapping(ZFLOW_LIQUID_CONDUCTANCE_INDEX)
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
              zflow_auxvar_up%pres - option%flow%reference_pressure < eps) then
            delta_pressure = 0.d0
            ddelta_pressure_dpdn = 0.d0
          endif
        endif
        if (delta_pressure >= 0.d0) then
          kr = zflow_auxvar_up%kr
          dkr_dp = 0.d0
        else
          kr = zflow_auxvar_dn%kr
          dkr_dp = zflow_auxvar_dn%dkr_dp
        endif

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(1) = perm_ave_over_dist_visc * kr * delta_pressure
      endif
    case(NEUMANN_BC)
      idof = auxvar_mapping(ZFLOW_LIQUID_FLUX_INDEX)
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(1) = auxvars(idof)
      endif
      derivative_toggle = PETSC_FALSE
    case default
      option%io_buffer = &
        'Boundary condition type (' // trim(StringWrite(bc_type)) // &
        ') not recognized in ZFlowBCFlux phase loop.'
      call PrintErrMsg(option)
  end select
  if (dabs(v_darcy(1)) > 0.d0 .or. kr > 0.d0) then
    ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(1) * area * zflow_density_kmol
    ! Res[m^3 liquid/sec]
    Res(1) = Res(1) + q
    if (calculate_derivatives .and. derivative_toggle) then
      ! derivative toggle takes care of the NEUMANN side
      tempreal = area * zflow_density_kmol * &
                 (dkr_dp * delta_pressure + ddelta_pressure_dpdn * kr)
      Jdn(1,1) = perm_ave_over_dist_visc * tempreal
      if (zflow_calc_adjoint) then
        if (zflow_adjoint_parameter == ZFLOW_ADJOINT_PERMEABILITY) then
          dJdndKdn(1) = dperm_dK * tempreal
          tempreal = area * zflow_density_kmol * kr
          drhsdKdn(1) = tempreal * (boundary_pressure + gravity_term) * dperm_dK
          dResdKdn(1) = Res(1) / perm_ave_over_dist_visc * dperm_dK
        endif
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
  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux(1)
  PetscReal :: scale
  PetscReal :: Res(1)
  PetscReal :: Jdn(1,1)
  PetscBool :: calculate_derivatives

  PetscReal :: qsrc_m3

  Res = 0.d0
  Jdn = 0.d0

  ! liquid phase
  qsrc_m3 = 0.d0
  select case(flow_src_sink_type)
    case(VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec
      qsrc_m3 = qsrc(1)
    case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
      qsrc_m3 = qsrc(1)*scale
    case default
      option%io_buffer = 'src_sink_type not supported in ZFlowSrcSink'
      call PrintErrMsg(option)
  end select
  ss_flow_vol_flux(1) = qsrc_m3
  ! Res[m^3 liquid/sec]
  Res = qsrc_m3 * zflow_density_kmol

  if (calculate_derivatives) then
  endif

end subroutine ZFlowSrcSink

! ************************************************************************** !

subroutine ZFlowAccumDerivative(zflow_auxvar,global_auxvar, &
                                material_auxvar, &
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

  type(zflow_auxvar_type) :: zflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: J(1,1)

  PetscReal :: res(1), res_pert(1)
  PetscReal :: Jdum(1,1)
  PetscReal :: dJdum(1,1)

  call ZFlowAccumulation(zflow_auxvar(ZERO_INTEGER), &
                         global_auxvar, &
                         material_auxvar, &
                         option,res,Jdum,dJdum, &
                         PETSC_FALSE)

  call ZFlowAccumulation(zflow_auxvar(ONE_INTEGER), &
                         global_auxvar, &
                         material_auxvar, &
                         option,res_pert,Jdum,dJdum, &
                         PETSC_FALSE)
  ! J[m^3 liquid/Pa-sec]
  J(1,1) = (res_pert(1)-res(1))/zflow_auxvar(ONE_INTEGER)%pert

end subroutine ZFlowAccumDerivative

! ************************************************************************** !

subroutine XXFluxDerivative(zflow_auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            zflow_auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            area, dist, &
                            zflow_parameter, &
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

  type(zflow_auxvar_type) :: zflow_auxvar_up(0:), zflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: Jup(1,1)
  PetscReal :: Jdn(1,1)

  PetscReal :: v_darcy(1)
  PetscReal :: res_up(1), res_dn(1)
  PetscReal :: res_pert(1)
  PetscReal :: Jdum(1,1), dJdum(1)

  Jup = 0.d0
  Jdn = 0.d0

  call XXFlux(zflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist, &
              zflow_parameter, &
              option,v_darcy,res_up,Jdum,Jdum, &
              dJdum,dJdum,dJdum,dJdum, &
              dJdum,dJdum, &
              dJdum,dJdum, &
              PETSC_FALSE)
  res_dn = res_up

  ! upgradient derivatives
  call XXFlux(zflow_auxvar_up(ONE_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist, &
              zflow_parameter, &
              option,v_darcy,res_pert,Jdum,Jdum, &
              dJdum,dJdum,dJdum,dJdum, &
              dJdum,dJdum, &
              dJdum,dJdum, &
              PETSC_FALSE)
  ! J[m^3 liquid/Pa-sec]
  Jup(1,1) = (res_pert(1)-res_up(1)) / &
             zflow_auxvar_up(ONE_INTEGER)%pert
  ! downgradient derivatives
  call XXFlux(zflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              zflow_auxvar_dn(ONE_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist, &
              zflow_parameter, &
              option,v_darcy,res_pert,Jdum,Jdum, &
              dJdum,dJdum,dJdum,dJdum, &
              dJdum,dJdum, &
              dJdum,dJdum, &
              PETSC_FALSE)
  ! J[m^3 liquid/Pa-sec]
  Jdn(1,1) = (res_pert(1)-res_dn(1)) / &
             zflow_auxvar_dn(ONE_INTEGER)%pert

end subroutine XXFluxDerivative

! ************************************************************************** !

subroutine XXBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                              zflow_auxvar_up, &
                              global_auxvar_up, &
                              zflow_auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              area,dist, &
                              zflow_parameter, &
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
  PetscInt :: ibndtype(1:1)
  PetscInt :: auxvar_mapping(ZFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(zflow_auxvar_type) :: zflow_auxvar_up, zflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(zflow_parameter_type) :: zflow_parameter
  PetscReal :: Jdn(1,1)

  PetscReal :: v_darcy(1)
  PetscReal :: res(1), res_pert(1)
  PetscReal :: Jdum(1,1), dJdum(1)

  Jdn = 0.d0

  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                zflow_auxvar_up,global_auxvar_up, &
                zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                material_auxvar_dn, &
                area,dist, &
                zflow_parameter, &
                option,v_darcy,res,Jdum,dJdum,dJdum,dJdum, &
                PETSC_FALSE)

  ! downgradient derivatives
  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                zflow_auxvar_up,global_auxvar_up, &
                zflow_auxvar_dn(ONE_INTEGER),global_auxvar_dn, &
                material_auxvar_dn, &
                area,dist, &
                zflow_parameter, &
                option,v_darcy,res_pert,Jdum,dJdum,dJdum,dJdum, &
                PETSC_FALSE)
  ! J[m^3 liquid/Pa-sec]
  Jdn(1,1) = (res_pert(1)-res(1))/zflow_auxvar_dn(ONE_INTEGER)%pert

end subroutine XXBCFluxDerivative

! ************************************************************************** !

subroutine ZFlowSrcSinkDerivative(option,qsrc,flow_src_sink_type, &
                                  zflow_auxvars,global_auxvar, &
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
  type(zflow_auxvar_type) :: zflow_auxvars(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: scale
  PetscReal :: Jac(1,1)

  PetscReal :: res(1), res_pert(1)
  PetscReal :: dummy_real(1)
  PetscReal :: Jdum(1,1)

  Jac = 0.d0
  ! unperturbed zflow_auxvars value
  call ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                    zflow_auxvars(ZERO_INTEGER),global_auxvar, &
                    material_auxvar,dummy_real,scale,res,Jdum, &
                    PETSC_FALSE)

  ! perturbed zflow_auxvars values
  call ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                    zflow_auxvars(ONE_INTEGER),global_auxvar, &
                    material_auxvar,dummy_real, &
                    scale,res_pert,Jdum, &
                    PETSC_FALSE)
  ! J[m^3 liquid/Pa-sec]
  Jac(1,1) = (res_pert(1)-res(1))/zflow_auxvars(ONE_INTEGER)%pert

end subroutine ZFlowSrcSinkDerivative

end module ZFlow_Common_module
