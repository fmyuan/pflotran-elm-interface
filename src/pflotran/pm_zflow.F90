module PM_ZFlow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module
  use ZFlow_Aux_module

  implicit none

  private

  PetscInt, parameter :: MAX_CHANGE_LIQ_PRES_NI = 1
  PetscInt, parameter :: MAX_CHANGE_CONC_NI = 2
  PetscInt, parameter :: MAX_RES_LIQ_EQ = 3
  PetscInt, parameter :: MAX_RES_SOL_EQ = 4

  type, public, extends(pm_subsurface_flow_type) :: pm_zflow_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: max_allow_liq_pres_change_ni
    PetscReal :: liq_pres_change_ts_governor
    PetscReal :: liq_sat_change_ts_governor
    PetscInt :: convergence_flags(MAX_RES_SOL_EQ)
    PetscReal :: convergence_reals(MAX_RES_SOL_EQ)
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMZFlowReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMZFlowReadTSSelectCase
    procedure, public :: ReadNewtonBlock => PMZFlowReadNewtonSelectCase
    procedure, public :: InitializeRun => PMZFlowInitializeRun
    procedure, public :: InitializeTimestep => PMZFlowInitializeTimestep
    procedure, public :: Residual => PMZFlowResidual
    procedure, public :: Jacobian => PMZFlowJacobian
    procedure, public :: UpdateTimestep => PMZFlowUpdateTimestep
    procedure, public :: FinalizeTimestep => PMZFlowFinalizeTimestep
    procedure, public :: PreSolve => PMZFlowPreSolve
    procedure, public :: PostSolve => PMZFlowPostSolve
    procedure, public :: CheckUpdatePre => PMZFlowCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMZFlowCheckUpdatePost
    procedure, public :: CheckConvergence => PMZFlowCheckConvergence
    procedure, public :: TimeCut => PMZFlowTimeCut
    procedure, public :: UpdateSolution => PMZFlowUpdateSolution
    procedure, public :: UpdateAuxVars => PMZFlowUpdateAuxVars
    procedure, public :: MaxChange => PMZFlowMaxChange
    procedure, public :: ComputeMassBalance => PMZFlowComputeMassBalance
    procedure, public :: InputRecord => PMZFlowInputRecord
    procedure, public :: CheckpointBinary => PMZFlowCheckpointBinary
    procedure, public :: RestartBinary => PMZFlowRestartBinary
    procedure, public :: Destroy => PMZFlowDestroy
  end type pm_zflow_type

  public :: PMZFlowCreate, &
            PMZFlowInitObject, &
            PMZFlowInitializeRun, &
            PMZFlowFinalizeTimestep, &
            PMZFlowCheckUpdatePre, &
            PMZFlowDestroy

contains

! ************************************************************************** !

function PMZFlowCreate()
  !
  ! Creates ZFlow process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  implicit none

  class(pm_zflow_type), pointer :: PMZFlowCreate

  class(pm_zflow_type), pointer :: zflow_pm

  allocate(zflow_pm)
  call PMZFlowInitObject(zflow_pm)

  PMZFlowCreate => zflow_pm

end function PMZFlowCreate

! ************************************************************************** !

subroutine PMZFlowInitObject(this)
  !
  ! Creates ZFlow process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Option_module
  use String_module
  use Variables_module

  implicit none

  class(pm_zflow_type) :: this

  call PMSubsurfaceFlowInit(this)
  this%name = 'Z Flow'
  this%header = 'Z FLOW'

  nullify(this%max_change_ivar)

  ! set to UNINITIALIZED_DOUBLE and report error below is set from input
  this%pressure_change_governor = UNINITIALIZED_DOUBLE
  this%temperature_change_governor = UNINITIALIZED_DOUBLE
  this%saturation_change_governor = UNINITIALIZED_DOUBLE
  this%xmol_change_governor = UNINITIALIZED_DOUBLE

  this%check_post_convergence = PETSC_TRUE

  this%max_allow_liq_pres_change_ni = 1.d20
  this%liq_pres_change_ts_governor = 5.d5    ! [Pa]
  this%liq_sat_change_ts_governor = 1.d0
  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMZFlowInitObject

! ************************************************************************** !

subroutine PMZFlowReadSimOptionsBlock(this,input)
  !
  ! Read ZFLOW options input block
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use ZFlow_module
  use ZFlow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module
  use EOS_Water_module

  implicit none

  class(pm_zflow_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH) :: local_error_string
  PetscBool :: found
  PetscReal :: array(1)

  option => this%option

  error_string = 'ZFlow Options'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('PROCESSES')
        local_error_string = trim(error_string) // ',' // keyword
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword',local_error_string)
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('LIQUID_FLOW')
              zflow_liq_flow_eq = 1     ! these serve as flags until updated
            case('HEAT_TRANSFER')
              zflow_heat_tran_eq = 1
            case('SOLUTE_TRANSPORT')
              zflow_sol_tran_eq = 1
            case default
              call InputKeywordUnrecognized(input,keyword,'ZFlow Mode',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('NO_ACCUMULATION')
        zflow_calc_accum = PETSC_FALSE
      case('NO_FLUX')
        zflow_calc_flux = PETSC_FALSE
      case('NO_BCFLUX')
        zflow_calc_bcflux = PETSC_FALSE
      case('TENSORIAL_RELATIVE_PERMEABILITY')
        zflow_tensorial_rel_perm = PETSC_TRUE
      case('LIQUID_DENSITY')
        call InputReadDouble(input,option,zflow_density_kg)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,zflow_density_kg,'kg/m^3', &
                                      trim(error_string)//','//keyword,option)
      case('LIQUID_VISCOSITY')
        call InputReadDouble(input,option,zflow_viscosity)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,zflow_viscosity,'Pa-s', &
                                      trim(error_string)//','//keyword,option)
      case default
        call InputKeywordUnrecognized(input,keyword,'ZFlow Mode',option)
    end select
  enddo
  call InputPopBlock(input,option)

  zflow_density_kmol = zflow_density_kg / FMWH2O
  array(1) = zflow_density_kg ! dist is the aux array
  call EOSWaterSetDensity('CONSTANT',array)

end subroutine PMZFlowReadSimOptionsBlock

! ************************************************************************** !

subroutine PMZFlowReadTSSelectCase(this,input,keyword,found, &
                                   error_string,option)
  !
  ! Read timestepper settings specific to the ZFLOW process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_zflow_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  found = PETSC_TRUE
  call PMSubsurfaceFlowReadTSSelectCase(this,input,keyword,found, &
                                        error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('LIQ_SAT_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_sat_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('LIQ_PRES_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_pres_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%liq_pres_change_ts_governor, &
                                    'Pa',keyword,option)
    case default
      found = PETSC_FALSE
  end select

end subroutine PMZFlowReadTSSelectCase

! ************************************************************************** !

subroutine PMZFlowReadNewtonSelectCase(this,input,keyword,found, &
                                       error_string,option)
  !
  ! Reads input file parameters associated with the ZFLOW process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use ZFlow_Aux_module

  implicit none

  class(pm_zflow_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: word

  error_string = 'ZFLOW Newton Solver'

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('REL_PERTURBATION')
      call InputReadDouble(input,option,zflow_rel_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      ! no units conversion since it is relative
    case('MIN_LIQ_PRESSURE_PERTURBATION')
      call InputReadDouble(input,option,zflow_pres_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      call InputReadAndConvertUnits(input,zflow_pres_min_pert, &
                                    'Pa',keyword,option)
    case('MAX_ALLOW_LIQ_PRES_CHANGE_NI')
      call InputReadDouble(input,option,this%max_allow_liq_pres_change_ni)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%max_allow_liq_pres_change_ni, &
                                    'Pa',keyword,option)
    case default
      found = PETSC_FALSE

  end select

end subroutine PMZFlowReadNewtonSelectCase

! ************************************************************************** !

recursive subroutine PMZFlowInitializeRun(this)
  !
  ! Initializes the ZFLOW mode run.
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Realization_Base_class
  use Patch_module
  use Field_module
  use Option_module
  use Variables_module

  implicit none

  class(pm_zflow_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscInt :: ivar

  patch => this%realization%patch
  field => this%realization%field
  option => this%option

  if (this%steady_state) zflow_calc_accum = PETSC_FALSE

  ivar = 0
  if (zflow_liq_flow_eq > 0) then
    ivar = ivar + 2
  endif
  if (zflow_sol_tran_eq > 0) then
    ivar = ivar + 1
  endif
  allocate(this%max_change_ivar(ivar))
  ivar = 1
  if (zflow_liq_flow_eq > 0) then
    this%max_change_ivar(ivar) = LIQUID_PRESSURE
    ivar = ivar + 1
    this%max_change_ivar(ivar) = LIQUID_SATURATION
    ivar = ivar + 1
  endif
  if (zflow_sol_tran_eq > 0) then
    this%max_change_ivar(ivar) = SOLUTE_CONCENTRATION
    ivar = ivar + 1
  endif

  ! need to allocate vectors for max change
  i = size(this%max_change_ivar)
  call VecDuplicateVecsF90(field%work,i,field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, size(field%max_change_vecs)
    call RealizationGetVariable(this%realization,field%max_change_vecs(i), &
                                this%max_change_ivar(i),ZERO_INTEGER)
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

  if (Initialized(this%temperature_change_governor)) then
    option%io_buffer = 'TEMPERATURE_CHANGE_GOVERNOR &
      &may not be used with ZFLOW.'
    call PrintErrMsg(option)
  endif

end subroutine PMZFlowInitializeRun

! ************************************************************************** !

subroutine PMZFlowInitializeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use ZFlow_module, only : ZFlowInitializeTimestep

  implicit none

  class(pm_zflow_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call ZFlowInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMZFlowInitializeTimestep

! ************************************************************************** !

subroutine PMZFlowFinalizeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  implicit none

  class(pm_zflow_type) :: this

  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMZFlowFinalizeTimestep

! ************************************************************************** !

subroutine PMZFlowPreSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  implicit none

  class(pm_zflow_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMZFlowPreSolve

! ************************************************************************** !

subroutine PMZFlowPostSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_zflow_type) :: this

end subroutine PMZFlowPostSolve

! ************************************************************************** !

subroutine PMZFlowUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                 num_newton_iterations,tfac, &
                                 time_step_max_growth_factor)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Option_module
  use Utility_module, only : Equal

  implicit none

  class(pm_zflow_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: sat_ratio, pres_ratio
  PetscReal :: dt_prev

  dt_prev = dt

  ! calculate the time step ramping factor
  sat_ratio = (2.d0*this%liq_sat_change_ts_governor)/ &
              (this%liq_sat_change_ts_governor+this%max_saturation_change)
  pres_ratio = (2.d0*this%liq_pres_change_ts_governor)/ &
               (this%liq_pres_change_ts_governor+this%max_pressure_change)
  ! pick minimum time step from calc'd ramping factor or maximum ramping factor
  dt = min(min(sat_ratio,pres_ratio)*dt,time_step_max_growth_factor*dt)
  ! make sure time step is within bounds given in the input deck
  dt = min(dt,dt_max)
  if (this%logging_verbosity > 0) then
    if (Equal(dt,dt_max)) then
      string = 'maximum time step size'
    else if (min(sat_ratio,pres_ratio) > time_step_max_growth_factor) then
      string = 'maximum time step growth factor'
    else if (sat_ratio < pres_ratio) then
      string = 'liquid saturation governor'
    else
      string = 'liquid pressure governor'
    endif
    string = 'TS update: ' // trim(string)
    call OptionPrint(string,this%option)
  endif

  if (Initialized(this%cfl_governor)) then
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMZFlowUpdateTimestep

! ************************************************************************** !

subroutine PMZFlowResidual(this,snes,xx,r,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use ZFlow_module, only : ZFlowResidual
  use Debug_module
  use Grid_module

  implicit none

  class(pm_zflow_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  PetscViewer :: viewer
  Mat :: M
  character(len=MAXSTRINGLENGTH) :: string

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  ! calculate residual
  if (zflow_simult_function_evals) then
    call SNESGetJacobian(snes,M,PETSC_NULL_MAT, &
                         PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER, &
                         ierr);CHKERRQ(ierr)
    call ZFlowResidual(snes,xx,r,M,this%realization,ierr)
  else
    call ZFlowResidual(snes,xx,r,PETSC_NULL_MAT,this%realization,ierr)
  endif

  if (this%realization%debug%vecview_residual) then
    string = 'ZFresidual'
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (this%realization%debug%vecview_solution) then
    string = 'ZFxx'
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine PMZFlowResidual

! ************************************************************************** !

subroutine PMZFlowJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Debug_module
  use Option_module

  implicit none

  class(pm_zflow_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: norm

  ! the Jacobian was already calculated in PMZFlowResidual

  if (this%realization%debug%matview_Matrix) then
    call DebugWriteFilename(this%realization%debug,string,'ZFjacobian','', &
                            zflow_ts_count,zflow_ts_cut_count, &
                            zflow_ni_count)
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  if (this%realization%debug%norm_Matrix) then
    call MatNorm(A,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("1 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("2 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("inf norm: ",es11.4)') norm
    call PrintMsg(this%option)
  endif

  zflow_ni_count = zflow_ni_count + 1

end subroutine PMZFlowJacobian

! ************************************************************************** !

subroutine PMZFlowCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  implicit none

  class(pm_zflow_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  changed = PETSC_FALSE

end subroutine PMZFlowCheckUpdatePre

! ************************************************************************** !

subroutine PMZFlowCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                  X1_changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_module
  use ZFlow_Aux_module

  implicit none

  class(pm_zflow_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscReal :: tempreal

  PetscBool :: converged_liquid_pressure
  PetscReal :: max_abs_pressure_change_NI
  PetscInt :: max_abs_pressure_change_NI_cell
  PetscBool :: converged_concentration
  PetscReal :: max_abs_conc_change_NI
  PetscInt :: max_abs_conc_change_NI_cell

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  ! If these are changed from true, we must add a global reduction on both
  ! variables to ensure that their values match across all processes. Otherwise
  ! PETSc will throw an error in debug mode or ignore the error in optimized.
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
  converged_liquid_pressure = PETSC_TRUE
  converged_concentration = PETSC_FALSE
  max_abs_pressure_change_NI = 0.d0
  max_abs_pressure_change_NI_cell = 0
  max_abs_conc_change_NI = 0.d0
  max_abs_conc_change_NI_cell = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    offset = (local_id-1)*option%nflowdof
    if (zflow_liq_flow_eq > 0) then
      ! maximum absolute change in liquid pressure over Newton iteration
      tempreal = dabs(dX_p(offset+zflow_liq_flow_eq))
      if (tempreal > dabs(max_abs_pressure_change_NI)) then
        max_abs_pressure_change_NI_cell = local_id
        max_abs_pressure_change_NI = tempreal
      endif
    endif
    if (zflow_sol_tran_eq > 0) then
      ! maximum absolute change in liquid pressure over Newton iteration
      tempreal = dabs(dX_p(offset+zflow_sol_tran_eq))
      if (tempreal > dabs(max_abs_conc_change_NI)) then
        max_abs_conc_change_NI_cell = local_id
        max_abs_conc_change_NI = tempreal
      endif
    endif
  enddo

  if (max_abs_pressure_change_NI > this%max_allow_liq_pres_change_ni) then
    converged_liquid_pressure = PETSC_FALSE
  endif

  ! the following flags are used in detemining convergence
  if (.not.converged_liquid_pressure) then
    this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = &
      max_abs_pressure_change_NI_cell
  endif

  ! the following flags are for REPORTING purposes only
  this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI) = max_abs_pressure_change_NI
  this%convergence_reals(MAX_CHANGE_CONC_NI) = max_abs_conc_change_NI

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)

end subroutine PMZFlowCheckUpdatePost

! ************************************************************************** !

subroutine PMZFlowCheckConvergence(this,snes,it,xnorm,unorm, &
                                   fnorm,reason,ierr)
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_module
  use ZFlow_Aux_module
  use Convergence_module

  implicit none

  class(pm_zflow_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  Vec :: residual_vec
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscReal, pointer :: X1_p(:)
  character(len=10) :: reason_string

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: converged_flag

  PetscReal, parameter :: zero_accumulation = 1.d-15

  PetscReal :: max_abs_res_liq_
  PetscInt :: max_abs_res_liq_cell
  PetscReal :: max_abs_res_sol_
  PetscInt :: max_abs_res_sol_cell
  PetscInt :: offset
  PetscMPIInt :: int_mpi

  PetscReal :: accumulation
  PetscReal :: residual
  PetscReal :: tempreal

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  zflow_auxvars => patch%aux%ZFlow%auxvars
  material_auxvars => patch%aux%Material%auxvars

  residual_vec = field%flow_r
  ! check residual terms
  call VecGetArrayReadF90(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (zflow_liq_flow_eq > 0) then
      residual = r_p(offset+zflow_liq_flow_eq)
      accumulation = accum2_p(offset+zflow_liq_flow_eq)
      ! residual
      tempreal = dabs(residual)
      if (tempreal > max_abs_res_liq_) then
        max_abs_res_liq_ = tempreal
        max_abs_res_liq_cell = local_id
      endif
    endif
    if (zflow_sol_tran_eq > 0) then
      residual = r_p(offset+zflow_sol_tran_eq)
      accumulation = accum2_p(offset+zflow_sol_tran_eq)
      ! residual
      tempreal = dabs(residual)
      if (tempreal > max_abs_res_sol_) then
        max_abs_res_sol_ = tempreal
        max_abs_res_sol_cell = local_id
      endif
    endif
  enddo

  ! the following flags are used in detemining convergence
  ! currently none

  ! the following flags are for REPORTING purposes only
  this%convergence_flags(MAX_RES_LIQ_EQ) = max_abs_res_liq_cell
  this%convergence_reals(MAX_RES_LIQ_EQ) = max_abs_res_liq_
  this%convergence_flags(MAX_RES_SOL_EQ) = max_abs_res_sol_cell
  this%convergence_reals(MAX_RES_SOL_EQ) = max_abs_res_sol_

  int_mpi = size(this%convergence_flags)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_flags,int_mpi, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  int_mpi = size(this%convergence_reals)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_reals,int_mpi, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

  ! these conditionals cannot change order
  reason_string = '---| '
  converged_flag = CONVERGENCE_CONVERGED
  if (this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) > 0) then
    reason_string(1:1) = 'P'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_CONC_NI) > 0) then
    reason_string(2:2) = 'C'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
#if 0
  if (OptionPrintToScreen(option)) then
    if (option%comm%mycommsize > 1 .or. grid%nmax > 9999) then
      write(*,'(4x,"Rsn: ",a10,2es10.2)') reason_string, &
        this%convergence_reals(MAX_RES_LIQ), &
        this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI)
    else
      write(*,'(4x,"Rsn: ",a10,2(i5,es10.2))') reason_string, &
        this%convergence_flags(MAX_RES_LIQ), &
        this%convergence_reals(MAX_RES_LIQ), &
        this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI), &
        this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI)
    endif
  endif
#endif
  option%convergence = converged_flag

  call VecRestoreArrayReadF90(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMZFlowCheckConvergence

! ************************************************************************** !

subroutine PMZFlowTimeCut(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use ZFlow_module, only : ZFlowTimeCut

  implicit none

  class(pm_zflow_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call ZFlowTimeCut(this%realization)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMZFlowTimeCut

! ************************************************************************** !

subroutine PMZFlowUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use ZFlow_module, only : ZFlowUpdateSolution, &
                           ZFlowMapBCAuxVarsToGlobal

  implicit none

  class(pm_zflow_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call ZFlowUpdateSolution(this%realization)
  call ZFlowMapBCAuxVarsToGlobal(this%realization)

end subroutine PMZFlowUpdateSolution

! ************************************************************************** !

subroutine PMZFlowUpdateAuxVars(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  use ZFlow_module, only : ZFlowUpdateAuxVars

  implicit none

  class(pm_zflow_type) :: this

  call ZFlowUpdateAuxVars(this%realization)

end subroutine PMZFlowUpdateAuxVars

! ************************************************************************** !

subroutine PMZFlowMaxChange(this)
  !
  ! Not needed given ZFlowMaxChange is called in PostSolve
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use String_module
  use ZFlow_Aux_module

  implicit none

  class(pm_zflow_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal, allocatable :: max_change_global(:)
  PetscReal :: max_change, change
  PetscInt :: i, j
  PetscInt :: ivar
  PetscInt :: fids(2)
  character(len=MAXSTRINGLENGTH) :: string

  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  allocate(max_change_global(size(this%max_change_ivar)))
  max_change_global = 0.d0

  do i = 1, size(this%max_change_ivar)
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i),ZERO_INTEGER)
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_old_ptr,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      change = dabs(vec_new_ptr(j)-vec_old_ptr(j))
      max_change = max(max_change,change)
    enddo
    max_change_global(i) = max_change
    call VecRestoreArrayF90(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_old_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  i = size(max_change_global)
  call MPI_Allreduce(MPI_IN_PLACE,max_change_global,i, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

  fids = OptionGetFIDs(option)
  ivar = 1
  if (zflow_liq_flow_eq > 0) then
    write(string,'("  --> max chng: dpl= ",1pe12.4, " dsl= ",1pe12.4)') &
      max_change_global(ivar:ivar+1)
    this%max_pressure_change = max_change_global(ivar)
    this%max_saturation_change = max_change_global(ivar+1)
    ivar = ivar+2
    call StringWriteToUnits(fids,string)
  endif
  if (zflow_sol_tran_eq > 0) then
    write(string,'("                 dc= ",1pe12.4)') max_change_global(ivar)
    ! hijacking xmol_change
    this%max_xmol_change = max_change_global(ivar)
    ivar = ivar+1
    call StringWriteToUnits(fids,string)
  endif
  deallocate(max_change_global)

end subroutine PMZFlowMaxChange

! ************************************************************************** !

subroutine PMZFlowComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use ZFlow_module, only : ZFlowComputeMassBalance

  implicit none

  class(pm_zflow_type) :: this
  PetscReal :: mass_balance_array(:)

  call ZFlowComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMZFlowComputeMassBalance

! ************************************************************************** !

subroutine PMZFlowInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 08/13/21
  !

  implicit none

  class(pm_zflow_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'zflow'

end subroutine PMZFlowInputRecord

! ************************************************************************** !

subroutine PMZFlowCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with ZFlow PM
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_zflow_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMZFlowCheckpointBinary

! ************************************************************************** !

subroutine PMZFlowRestartBinary(this,viewer)
  !
  ! Restarts data associated with ZFlow PM
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_zflow_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMZFlowRestartBinary

! ************************************************************************** !

subroutine PMZFlowDestroy(this)
  !
  ! Destroys ZFlow process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use ZFlow_module, only : ZFlowDestroy
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_zflow_type) :: this

  PetscErrorCode :: ierr

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call DeallocateArray(this%max_change_ivar)
  call ZFlowDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMZFlowDestroy

end module PM_ZFlow_class
