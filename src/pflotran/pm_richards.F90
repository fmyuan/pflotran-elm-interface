module PM_Richards_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: ABS_UPDATE_INDEX = 1
  PetscInt, parameter :: REL_UPDATE_INDEX = 2
  PetscInt, parameter :: RESIDUAL_INDEX = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  PetscInt, parameter :: MAX_INDEX = SCALED_RESIDUAL_INDEX

  type, public, extends(pm_subsurface_flow_type) :: pm_richards_type
    PetscBool :: converged_flag(MAX_INDEX)
    PetscInt :: converged_cell(MAX_INDEX)
    PetscReal :: converged_real(MAX_INDEX)
    PetscReal :: residual_abs_inf_tol
    PetscReal :: residual_scaled_inf_tol
    PetscReal :: abs_update_inf_tol
    PetscReal :: rel_update_inf_tol
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMRichardsReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock => PMRichardsReadNewtonSelectCase
    procedure, public :: InitializeTimestep => PMRichardsInitializeTimestep
    procedure, public :: Residual => PMRichardsResidual
    procedure, public :: Jacobian => PMRichardsJacobian
    procedure, public :: UpdateTimestep => PMRichardsUpdateTimestep
    procedure, public :: PreSolve => PMRichardsPreSolve
    procedure, public :: PostSolve => PMRichardsPostSolve
    procedure, public :: CheckUpdatePre => PMRichardsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRichardsCheckUpdatePost
    procedure, public :: CheckConvergence => PMRichardsCheckConvergence
    procedure, public :: TimeCut => PMRichardsTimeCut
    procedure, public :: UpdateSolution => PMRichardsUpdateSolution
    procedure, public :: UpdateAuxVars => PMRichardsUpdateAuxVars
    procedure, public :: MaxChange => PMRichardsMaxChange
    procedure, public :: ComputeMassBalance => PMRichardsComputeMassBalance
    procedure, public :: InputRecord => PMRichardsInputRecord
    procedure, public :: Destroy => PMRichardsDestroy
  end type pm_richards_type

  public :: PMRichardsCreate, &
            PMRichardsInit, &
            PMRichardsDestroy, &
            PMRichardsCheckConvergence


contains

! ************************************************************************** !

function PMRichardsCreate()
  !
  ! Creates Richards process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  implicit none

  class(pm_richards_type), pointer :: PMRichardsCreate

  class(pm_richards_type), pointer :: this

  allocate(this)
  call PMRichardsInit(this)
  this%name = 'Richards Flow'
  this%header = 'RICHARDS FLOW'

  this%residual_abs_inf_tol = 1.d-5
  this%residual_scaled_inf_tol = 1.d-5
  this%abs_update_inf_tol = 1.d0
  this%rel_update_inf_tol = 1.d-5

  PMRichardsCreate => this

end function PMRichardsCreate

! ************************************************************************** !

subroutine PMRichardsInit(this)
  !
  ! Initializes Richards process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 07/11/18
  !

  implicit none

  class(pm_richards_type) :: this

  call PMSubsurfaceFlowInit(this)

end subroutine PMRichardsInit

! ************************************************************************** !

subroutine PMRichardsReadSimOptionsBlock(this,input)
  !
  ! Reads input file parameters associated with the Richards process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/16/20

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Richards_Aux_module

  implicit none

  class(pm_richards_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found
  PetscReal :: tempreal

  option => this%option

  error_string = 'Richards Options'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('INLINE_SURFACE_REGION')
        option%flow%inline_surface_flow = PETSC_TRUE
        call InputReadWord(input,option,keyword,PETSC_FALSE)
        option%flow%inline_surface_region_name = keyword
      case('INLINE_SURFACE_MANNINGS_COEFF')
        call InputReadDouble(input,option,tempreal)
        option%flow%inline_surface_Mannings_coeff = tempreal
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMRichardsReadSimOptionsBlock

! ************************************************************************** !

subroutine PMRichardsReadNewtonSelectCase(this,input,keyword,found, &
                                          error_string,option)
  !
  ! Reads input file parameters associated with the Richards process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 03/16/20

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Richards_Aux_module

  implicit none

  class(pm_richards_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  PetscBool :: found
  PetscReal :: tempreal

  option => this%option

  error_string = 'Richards Newton Solver'

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
! Tolerances
    ! All Residual
    case('RESIDUAL_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_abs_inf_tol = tempreal
      this%residual_scaled_inf_tol = tempreal
    ! Absolute Residual
    case('RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_abs_inf_tol = tempreal
    ! Scaled Residual
    case('RESIDUAL_SCALED_INF_TOL','ITOL_SCALED_RESIDUAL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_scaled_inf_tol = tempreal
    ! All Updates
    case('UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol = tempreal
      this%rel_update_inf_tol = tempreal
    ! Absolute Updates
    case('ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol = tempreal
    ! Relative Updates
    case('REL_UPDATE_INF_TOL','ITOL_RELATIVE_UPDATE')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol = tempreal
    case default
      found = PETSC_FALSE

  end select

end subroutine PMRichardsReadNewtonSelectCase

! ************************************************************************** !

subroutine PMRichardsInitializeTimestep(this)
  !
  ! Should not need this as it is called in PreSolve.
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsInitializeTimestep
  use Option_module

  implicit none

  class(pm_richards_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  call RichardsInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMRichardsInitializeTimestep

! ************************************************************************** !

subroutine PMRichardsPreSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_richards_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMRichardsPreSolve

! ************************************************************************** !

subroutine PMRichardsPostSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_richards_type) :: this

end subroutine PMRichardsPostSolve

! ************************************************************************** !

subroutine PMRichardsUpdateTimestep(this,update_dt, &
                                    dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac, &
                                    time_step_max_growth_factor)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none

  class(pm_richards_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac

  if (update_dt .and. iacceleration /= 0) then
    if (iacceleration > 0) then
      fac = 0.5d0
      if (num_newton_iterations >= iacceleration) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = this%pressure_change_governor/(this%max_pressure_change+0.1)
        ut = up
      endif
      dtt = fac * dt * (1.d0 + ut)
    else
      ifac = max(min(num_newton_iterations,size(tfac)),1)
      dt_tfac = tfac(ifac) * dt

      fac = 0.5d0
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      dt_p = fac * dt * (1.d0 + up)

      dtt = min(dt_tfac,dt_p)
    endif

    dtt = min(time_step_max_growth_factor*dt,dtt)
    if (dtt > dt_max) dtt = dt_max
    ! geh: There used to be code here that cut the time step if it is too
    !      large relative to the simulation time.  This has been removed.
    dtt = max(dtt,dt_min)
    dt = dtt
  endif

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)

end subroutine PMRichardsUpdateTimestep

! ************************************************************************** !

subroutine PMRichardsResidual(this,snes,xx,r,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsResidual

  implicit none

  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call RichardsResidual(snes,xx,r,this%realization,ierr)

end subroutine PMRichardsResidual

! ************************************************************************** !

subroutine PMRichardsJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsJacobian

  implicit none

  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call RichardsJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMRichardsJacobian

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Characteristic_Curves_module
  use Patch_module
  use Richards_Aux_module
  use Global_Aux_module
  use Patch_module

  implicit none

  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscReal :: P_R, P0, P1, delP
  PetscReal :: scale, sat, sat_pert, pert, pc_pert, press_pert, delP_pert
  PetscReal :: dpc_dsatl

  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars

  if (Initialized(this%saturation_change_limit)) then

    changed = PETSC_TRUE

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)

    pert = dabs(this%saturation_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      sat = global_auxvars(ghosted_id)%sat(1)
      sat_pert = sat - sign(1.d0,sat-0.5d0)*pert
      call patch%characteristic_curves_array( &
        patch%cc_id(ghosted_id))%ptr%saturation_function% &
        CapillaryPressure(sat_pert,pc_pert,dpc_dsatl,option)
      press_pert = option%flow%reference_pressure - pc_pert
      P0 = X_p(local_id)
      delP = dX_p(local_id)
      delP_pert = dabs(P0 - press_pert)
      if (delP_pert < dabs(delP)) then
        write(option%io_buffer,'("dP_trunc:",1i7,2es15.7)') &
          grid%nG2A(grid%nL2G(local_id)),delP_pert,dabs(delP)
        call PrintMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),delP_pert),delP)
      dX_p(local_id) = delP
    enddo

    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)

  endif

  if (Initialized(this%pressure_dampening_factor)) then
    changed = PETSC_TRUE
    ! P^p+1 = P^p - dP^p
    P_R = option%flow%reference_pressure
    scale = this%pressure_dampening_factor

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      delP = dX_p(local_id)
      P0 = X_p(local_id)
      P1 = P0 - delP
      if (P0 < P_R .and. P1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call PrintMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call PrintMsgAnyRank(option)
#endif
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call PrintMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call PrintMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dX_p(local_id) = scale*delP
      endif
    enddo
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine PMRichardsCheckUpdatePre

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                     X1_changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/21/18
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none

  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: dX_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  PetscInt :: local_id, ghosted_id, natural_id, ival
  PetscReal :: dX_, X0_, dX_X0
  PetscBool :: converged_abs_update_flag
  PetscBool :: converged_rel_update_flag
  PetscInt :: converged_abs_update_cell
  PetscInt :: converged_rel_update_cell
  PetscReal :: converged_abs_update_real
  PetscReal :: converged_rel_update_real
  PetscBool :: converged_absolute
  PetscBool :: converged_relative

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  converged_abs_update_flag = PETSC_TRUE
  converged_rel_update_flag = PETSC_TRUE
  converged_abs_update_cell = ZERO_INTEGER
  converged_rel_update_cell = ZERO_INTEGER
  converged_abs_update_real = 0.d0
  converged_rel_update_real = 0.d0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    ival = local_id
    ! infinity norms on update
    converged_absolute = PETSC_TRUE
    converged_relative = PETSC_TRUE
    dX_ = dabs(dX_p(ival))
    X0_ = dabs(X0_p(ival))
    X0_ = max(X0_,1.d-40)
    dX_X0 = dX_/X0_
    if (dX_ > this%abs_update_inf_tol) then
      converged_absolute = PETSC_FALSE
    endif
    if (converged_abs_update_real < dX_) then
      converged_abs_update_real = dX_
      converged_abs_update_cell = natural_id
    endif
    if (dX_X0 > this%rel_update_inf_tol) then
      converged_relative = PETSC_FALSE
    endif
    if (converged_rel_update_real < dX_X0) then
      converged_rel_update_real = dX_X0
      converged_rel_update_cell = natural_id
    endif
    ! only enter this condition if both are not converged
    if (.not.(converged_absolute .or. converged_relative)) then
      if (.not.converged_absolute) then
        converged_abs_update_flag = PETSC_FALSE
      endif
      if (.not.converged_relative) then
        converged_rel_update_flag = PETSC_FALSE
      endif
    endif
  enddo
  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(ABS_UPDATE_INDEX) = converged_abs_update_flag
  this%converged_flag(REL_UPDATE_INDEX) = converged_rel_update_flag
  this%converged_real(ABS_UPDATE_INDEX) = converged_abs_update_real
  this%converged_real(REL_UPDATE_INDEX) = converged_rel_update_real
  this%converged_cell(ABS_UPDATE_INDEX) = converged_abs_update_cell
  this%converged_cell(REL_UPDATE_INDEX) = converged_rel_update_cell

end subroutine PMRichardsCheckUpdatePost

! ************************************************************************** !

subroutine PMRichardsCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                      reason,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/20/18
  !
  use Convergence_module
  use General_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use String_module

  implicit none

  class(pm_richards_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: itol
  PetscReal :: R, A, R_A
  PetscReal, parameter :: A_zero = 1.d-15
  PetscBool :: converged_abs_residual_flag
  PetscReal :: converged_abs_residual_real
  PetscInt :: converged_abs_residual_cell
  PetscBool :: converged_scaled_residual_flag
  PetscReal :: converged_scaled_residual_real
  PetscInt :: converged_scaled_residual_cell
  PetscBool :: converged_absolute
  PetscBool :: converged_scaled
  PetscMPIInt :: mpi_int
  character(len=MAXSTRINGLENGTH) :: string
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid

  if (this%check_post_convergence) then
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
    converged_abs_residual_flag = PETSC_TRUE
    converged_abs_residual_real = 0.d0
    converged_abs_residual_cell = ZERO_INTEGER
    converged_scaled_residual_flag = PETSC_TRUE
    converged_scaled_residual_real = 0.d0
    converged_scaled_residual_cell = ZERO_INTEGER
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      converged_absolute = PETSC_TRUE
      converged_scaled = PETSC_TRUE
      ! infinity norms on residual
      R = dabs(r_p(local_id))
      A = dabs(accum2_p(local_id))
      R_A = R/A
      if (R > this%residual_abs_inf_tol) then
        converged_absolute = PETSC_FALSE
      endif
      ! find max value regardless of convergence
      if (converged_abs_residual_real < R) then
        converged_abs_residual_real = R
        converged_abs_residual_cell = natural_id
      endif
      if (A > A_zero) then
        if (R_A > this%residual_scaled_inf_tol) then
          converged_scaled = PETSC_FALSE
        endif
        ! find max value regardless of convergence
        if (converged_scaled_residual_real < R_A) then
          converged_scaled_residual_real = R_A
          converged_scaled_residual_cell = natural_id
        endif
      endif
      ! only enter this condition if both are not converged
      if (.not.(converged_absolute .or. converged_scaled)) then
        if (.not.converged_absolute) then
          converged_abs_residual_flag = PETSC_FALSE
        endif
        if (.not.converged_scaled) then
          converged_scaled_residual_flag = PETSC_FALSE
        endif
      endif
    enddo
    call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum2,accum2_p, &
                                ierr);CHKERRQ(ierr)

    this%converged_flag(RESIDUAL_INDEX) = converged_abs_residual_flag
    this%converged_real(RESIDUAL_INDEX) = converged_abs_residual_real
    this%converged_cell(RESIDUAL_INDEX) = converged_abs_residual_cell
    this%converged_flag(SCALED_RESIDUAL_INDEX) = converged_scaled_residual_flag
    this%converged_real(SCALED_RESIDUAL_INDEX) = converged_scaled_residual_real
    this%converged_cell(SCALED_RESIDUAL_INDEX) = converged_scaled_residual_cell
    mpi_int = MAX_INDEX
    ! do not perform an all reduce on cell id as this info is not printed
    ! in parallel
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_flag,mpi_int,MPI_INTEGER, &
                       MPI_LAND,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                       ierr);CHKERRQ(ierr)

    option%convergence = CONVERGENCE_CONVERGED

    do itol = 1, MAX_INDEX
      if (.not.this%converged_flag(itol)) then
        option%convergence = CONVERGENCE_KEEP_ITERATING
        if (this%logging_verbosity > 0) then
          string = '   ' // trim(tol_string(itol)) // ', Liquid Pressure'
          if (option%comm%mycommsize == 1) then
            string = trim(string) // ' (' // &
              trim(StringFormatInt(this%converged_cell(itol))) &
              // ')'
          endif
          string = trim(string) // ' : ' // &
            StringFormatDouble(this%converged_real(itol))
          call PrintMsg(option,string)
        endif
      endif
    enddo
    if (this%logging_verbosity > 0 .and. it > 0 .and. &
        option%convergence == CONVERGENCE_CONVERGED) then
      string = '   Converged'
      call PrintMsg(option,string)
      write(string,'(4x," R:",es8.1)') this%converged_real(RESIDUAL_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"SR:",es8.1)') &
        this%converged_real(SCALED_RESIDUAL_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"AU:",es8.1)') this%converged_real(ABS_UPDATE_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"RU:",es8.1)') this%converged_real(REL_UPDATE_INDEX)
      call PrintMsg(option,string)
    endif
  endif

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMRichardsCheckConvergence

! ************************************************************************** !

subroutine PMRichardsTimeCut(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsTimeCut

  implicit none

  class(pm_richards_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call RichardsTimeCut(this%realization)
  call PMSubsurfaceFlowTimeCutPostInit(this)

end subroutine PMRichardsTimeCut

! ************************************************************************** !

subroutine PMRichardsUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsUpdateSolution

  implicit none

  class(pm_richards_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call RichardsUpdateSolution(this%realization)

end subroutine PMRichardsUpdateSolution

! ************************************************************************** !

subroutine PMRichardsUpdateAuxVars(this)
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Richards_module, only : RichardsUpdateAuxVars

  implicit none

  class(pm_richards_type) :: this

  call RichardsUpdateAuxVars(this%realization)

end subroutine PMRichardsUpdateAuxVars

! ************************************************************************** !

subroutine PMRichardsMaxChange(this)
  !
  ! Not needed given RichardsMaxChange is called in PostSolve
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsMaxChange
  use Option_module

  implicit none

  class(pm_richards_type) :: this

  call RichardsMaxChange(this%realization,this%max_pressure_change)
  write(this%option%io_buffer,'("  --> max change: dpmx= ",1pe12.4)') &
        this%max_pressure_change
  call PrintMsg(this%option)

end subroutine PMRichardsMaxChange

! ************************************************************************** !

subroutine PMRichardsComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsComputeMassBalance

  implicit none

  class(pm_richards_type) :: this
  PetscReal :: mass_balance_array(:)

  call RichardsComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMRichardsComputeMassBalance

! ************************************************************************** !

subroutine PMRichardsInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  !

  implicit none

  class(pm_richards_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'richards'
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_UPDATE: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMRichardsInputRecord

! ************************************************************************** !

subroutine PMRichardsDestroy(this)
  !
  ! Destroys Richards process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Richards_module, only : RichardsDestroy

  implicit none

  class(pm_richards_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call RichardsDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMRichardsDestroy

end module PM_Richards_class
