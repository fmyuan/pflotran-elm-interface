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
                         option,v_darcy,Res, &
                         derivative_call, &
                         debug_connection)
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
      PetscBool :: derivative_call
      PetscBool :: debug_connection
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                         zflow_auxvar_up,global_auxvar_up, &
                         zflow_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist, &
                         zflow_parameter, &
                         option,v_darcy,Res, &
                         derivative_call, &
                         debug_connection)
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
      PetscBool :: derivative_call
      PetscBool :: debug_connection
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
                             option,Res,debug_cell)
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
  PetscBool :: debug_cell

  PetscReal :: porosity
  PetscReal :: volume_over_dt

  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use zflow_auxvar%effective porosity here as it enables numerical
  ! derivatives to be employed
  porosity = zflow_auxvar%effective_porosity

  ! accumulation term units = m^3 liquid/s
  ! Res[m^3 liquid/sec] = sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
  !                       vol[m^3 bulk] / dt[sec]
  Res(1) = zflow_auxvar%sat * porosity * volume_over_dt

end subroutine ZFlowAccumulation

! ************************************************************************** !

subroutine ZFlowFluxHarmonicPermOnly(zflow_auxvar_up,global_auxvar_up, &
                                   material_auxvar_up, &
                                   zflow_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   area, dist, &
                                   zflow_parameter, &
                                   option,v_darcy,Res, &
                                   derivative_call, &
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
  PetscBool :: derivative_call
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscReal :: perm_ave_over_dist
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q

  Res = 0.d0
  v_darcy = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  perm_ave_over_dist = (perm_up * perm_dn) / &
                       (dist_up*perm_dn + dist_dn*perm_up)

  if (zflow_auxvar_up%kr + zflow_auxvar_dn%kr > eps) then

    gravity_term = zflow_density_kg * dist_gravity
    delta_pressure = zflow_auxvar_up%pres - &
                     zflow_auxvar_dn%pres + &
                     gravity_term
    if (delta_pressure >= 0.d0) then
      kr = zflow_auxvar_up%kr
    else
      kr = zflow_auxvar_dn%kr
    endif

    if (kr > floweps ) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(1) = perm_ave_over_dist * kr * delta_pressure
      ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(1) * area
      ! Res[m^3 liquid/sec]
      Res = Res + q
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
                                     option,v_darcy,Res, &
                                     derivative_call, &
                                     debug_connection)
  !
  ! Computes the boundary flux terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  !
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
  PetscBool :: derivative_call
  PetscBool :: debug_connection

  PetscInt :: bc_type
  PetscInt :: idof
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: perm_dn
  PetscReal :: boundary_pressure
  PetscReal :: tempreal

  Res = 0.d0
  v_darcy = 0.d0

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  kr = 0.d0
  bc_type = ibndtype(1)
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
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn / dist(0)
        endif

        boundary_pressure = zflow_auxvar_up%pres
        gravity_term = zflow_density_kg * dist_gravity
        delta_pressure = boundary_pressure - &
                         zflow_auxvar_dn%pres + &
                         gravity_term
        if (bc_type == HYDROSTATIC_SEEPAGE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (delta_pressure > 0.d0 .and. &
              zflow_auxvar_up%pres - option%flow%reference_pressure < eps) then
            delta_pressure = 0.d0
          endif
        endif
        if (delta_pressure >= 0.d0) then
          kr = zflow_auxvar_up%kr
        else
          kr = zflow_auxvar_dn%kr
        endif

        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(1) = perm_ave_over_dist * kr * delta_pressure
      endif
    case(NEUMANN_BC)
      idof = auxvar_mapping(ZFLOW_LIQUID_FLUX_INDEX)
      !geh: we should read in the mole fraction for both phases as the
      !     enthalpy, etc. applies to phase, not pure component.
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy(1) = auxvars(idof)
      endif
    case default
      option%io_buffer = &
        'Boundary condition type not recognized in ZFlowBCFlux phase loop.'
      call PrintErrMsg(option)
  end select
  if (dabs(v_darcy(1)) > 0.d0 .or. kr > 0.d0) then
    ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
    q = v_darcy(1) * area
    ! Res[m^3 liquid/sec]
    Res(1) = Res(1) + q
  endif

end subroutine ZFlowBCFluxHarmonicPermOnly

! ************************************************************************** !

subroutine ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                        zflow_auxvar,global_auxvar,material_auxvar, &
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
  type(zflow_auxvar_type) :: zflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux(1)
  PetscReal :: scale
  PetscReal :: Res(1)
  PetscBool :: debug_cell

  PetscReal :: qsrc_m3
  PetscReal :: cell_pressure, dummy_pressure
  PetscErrorCode :: ierr

  Res = 0.d0

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
  Res = qsrc_m3

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
  PetscInt :: idof, irow

  call ZFlowAccumulation(zflow_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar,option, &
                           res,PETSC_FALSE)

  call ZFlowAccumulation(zflow_auxvar(ONE_INTEGER), &
                           global_auxvar, &
                           material_auxvar, &
                           option,res_pert,PETSC_FALSE)
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
  PetscInt :: idof

  Jup = 0.d0
  Jdn = 0.d0

  option%iflag = -2
  call XXFlux(zflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist, &
              zflow_parameter, &
              option,v_darcy,res_up, &
              PETSC_TRUE, & ! derivative call
              PETSC_FALSE)
  res_dn = res_up

  if (zflow_jacobian_test) then
    print *, 'res_dn: ', res_dn
    print *, 'res_up: ', res_up
  endif

  ! upgradient derivatives
  call XXFlux(zflow_auxvar_up(ONE_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist, &
              zflow_parameter, &
              option,v_darcy,res_pert, &
              PETSC_TRUE, & ! derivative call
              PETSC_FALSE)
  if (zflow_jacobian_test) then
    if (zflow_jacobian_test_xdof > 0) then
      print *, 'res_pert_up: ', res_pert
    endif
  endif
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
              option,v_darcy,res_pert, &
              PETSC_TRUE, & ! derivative call
              PETSC_FALSE)
  if (zflow_jacobian_test) then
    if (zflow_jacobian_test_xdof > 0) then
      print *, 'res_pert_dn: ', res_pert
    endif
  endif
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

  Jdn = 0.d0

  option%iflag = -2
  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     zflow_auxvar_up,global_auxvar_up, &
                     zflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist, &
                     zflow_parameter, &
                     option,v_darcy,res, &
                     PETSC_TRUE, & ! derivative call
                     PETSC_FALSE)

  ! downgradient derivatives
  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                      zflow_auxvar_up,global_auxvar_up, &
                      zflow_auxvar_dn(ONE_INTEGER),global_auxvar_dn, &
                      material_auxvar_dn, &
                      area,dist, &
                      zflow_parameter, &
                      option,v_darcy,res_pert, &
                      PETSC_TRUE, & ! derivative call
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

  option%iflag = -3
  ! unperturbed zflow_auxvars value
  call ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                      zflow_auxvars(ZERO_INTEGER),global_auxvar, &
                      material_auxvar,dummy_real,scale,res,PETSC_FALSE)

  ! perturbed zflow_auxvars values
  call ZFlowSrcSink(option,qsrc,flow_src_sink_type, &
                      zflow_auxvars(ONE_INTEGER),global_auxvar, &
                      material_auxvar,dummy_real, &
                      scale,res_pert,PETSC_FALSE)
  ! J[m^3 liquid/Pa-sec]
  Jac(1,1) = (res_pert(1)-res(1))/zflow_auxvars(ONE_INTEGER)%pert

end subroutine ZFlowSrcSinkDerivative

end module ZFlow_Common_module
