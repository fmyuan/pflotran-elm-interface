module Timestepper_Steady_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Solver_module
  use Convergence_module
  use Timestepper_SNES_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(timestepper_SNES_type) :: timestepper_steady_type

  contains

    procedure, public :: ReadSelectCase => TimestepperSteadyReadSelectCase
    procedure, public :: SetTargetTime => TimestepperSteadySetTargetTime
    procedure, public :: UpdateDT => TimestepperSteadyUpdateDT
    procedure, public :: StepDT => TimestepperSteadyStepDT
    procedure, public :: PrintInfo => TimestepperSteadyPrintInfo
    procedure, public :: InputRecord => TimestepperSteadyInputRecord
    procedure, public :: FinalizeRun => TimestepperSteadyFinalizeRun
    procedure, public :: Destroy => TimestepperSteadyDestroy

  end type timestepper_steady_type

  public :: TimestepperSteadyCreate, &
            TimestepperSteadyCast

contains

! ************************************************************************** !

function TimestepperSteadyCreate()
  !
  ! Allocates a new steady state timestepper object
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !

  implicit none

  class(timestepper_steady_type), pointer :: TimestepperSteadyCreate

  class(timestepper_steady_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  TimestepperSteadyCreate => stepper

end function TimestepperSteadyCreate

! ************************************************************************** !

function TimestepperSteadyCast(this)
  !
  ! Casts a timestepper object to its appropriate class
  !
  ! Author: Glenn Hammond
  ! Date: 03/29/21
  !
  use Timestepper_Base_class

  implicit none

  class(timestepper_base_type), pointer :: this

  class(timestepper_steady_type), pointer :: TimestepperSteadyCast

  nullify(TimestepperSteadyCast)
  if (.not.associated(this)) return
  select type(this)
    class is(timestepper_steady_type)
      TimestepperSteadyCast => this
  end select

end function TimestepperSteadyCast

! ************************************************************************** !

subroutine TimestepperSteadyReadSelectCase(this,input,keyword,found, &
                                           error_string,option)
  !
  ! Reads select case statement for Steady
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !
  use Option_module
  use Input_Aux_module

  implicit none

  class(timestepper_steady_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  option%io_buffer = 'TIMESTEPPER card not supported for steady state &
    &simulations. Please remove from ' // trim(error_string) // ' block.'
  call PrintErrMsg(option)

end subroutine TimestepperSteadyReadSelectCase

! ************************************************************************** !

subroutine TimestepperSteadyUpdateDT(this,process_model)
  !
  ! Updates time step
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !
  use PM_Base_class

  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model

  ! do nothing

end subroutine TimestepperSteadyUpdateDT

! ************************************************************************** !

subroutine TimestepperSteadySetTargetTime(this,sync_time,option,stop_flag, &
                                          sync_flag, &
                                          snapshot_plot_flag, &
                                          observation_plot_flag, &
                                          massbal_plot_flag,checkpoint_flag)
  !
  ! Sets target time for steady state ts, which
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !
  use Timestepper_Base_class
  use Option_module
  use Utility_module

  implicit none

  class(timestepper_steady_type) :: this
  PetscReal :: sync_time
  type(option_type) :: option
  PetscInt :: stop_flag
  PetscBool :: sync_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_flag

  this%target_time = sync_time
  do
    if (.not.associated(this%cur_waypoint)) exit
    if (sync_time >= this%cur_waypoint%time) then
      if (Equal(sync_time,this%cur_waypoint%time)) then
        if (this%cur_waypoint%print_snap_output) &
          sync_flag = PETSC_TRUE
        if (this%cur_waypoint%print_snap_output) &
          snapshot_plot_flag = PETSC_TRUE
        if (this%cur_waypoint%print_obs_output) &
          observation_plot_flag = PETSC_TRUE
        if (this%cur_waypoint%print_msbl_output) &
          massbal_plot_flag = PETSC_TRUE
        if (this%cur_waypoint%print_checkpoint) &
          checkpoint_flag = PETSC_TRUE
      endif
      this%cur_waypoint => this%cur_waypoint%next
    else
      exit
    endif
  enddo
  if (.not.associated(this%cur_waypoint)) then
    stop_flag = TS_STOP_END_SIMULATION
  endif

end subroutine TimestepperSteadySetTargetTime

! ************************************************************************** !

subroutine TimestepperSteadyStepDT(this,process_model,stop_flag)
  !
  ! Steps forward one step in time
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use Option_module
  use Output_module, only : Output, OutputFindNaNOrInfInVec
  use Output_EKG_module, only : IUNIT_EKG
  use String_module
  use Timestepper_Base_class

  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  SNESConvergedReason :: snes_reason

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option

  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations

  PetscReal :: fnorm, inorm, scaled_fnorm
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  Vec :: residual_vec
  PetscErrorCode :: ierr

  solver => process_model%solver
  option => process_model%option

  num_linear_iterations = 0
  num_newton_iterations = 0

  option%time = this%target_time

  call process_model%InitializeTimestep()

  call process_model%PreSolve()

  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

  call SNESSolve(solver%snes,PETSC_NULL_VEC,process_model%solution_vec, &
                 ierr);CHKERRQ(ierr)

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

  this%cumulative_solver_time = &
    this%cumulative_solver_time + &
    (log_end_time - log_start_time)

  call SNESGetIterationNumber(solver%snes,num_newton_iterations, &
                              ierr);CHKERRQ(ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, &
                                    ierr);CHKERRQ(ierr)
  call SNESGetConvergedReason(solver%snes,snes_reason,ierr);CHKERRQ(ierr)

  this%steps = this%steps + 1
  this%cumulative_newton_iterations = this%cumulative_newton_iterations + &
    num_newton_iterations
  this%cumulative_linear_iterations = this%cumulative_linear_iterations + &
    num_linear_iterations

  if (snes_reason <= SNES_CONVERGED_ITERATING .or. &
      .not. process_model%AcceptSolution()) then
    call SolverNewtonPrintFailedReason(solver,option)
    if (solver%verbose_logging) then
      select case(snes_reason)
        case(SNES_DIVERGED_FNORM_NAN)
          ! attempt to find cells with NaNs.
          call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_FUNCTION, &
                               PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
          call OutputFindNaNOrInfInVec(residual_vec, &
                                       process_model%realization_base% &
                                         discretization%grid,option)
      end select
    endif
    process_model%output_option%plot_name = trim(process_model%name) // &
      '_cut_to_failure'
    snapshot_plot_flag = PETSC_TRUE
    observation_plot_flag = PETSC_FALSE
    massbal_plot_flag = PETSC_FALSE
    call Output(process_model%realization_base,snapshot_plot_flag, &
                observation_plot_flag,massbal_plot_flag)
    option%io_buffer = 'Newton solver failed to converge in steady-state &
                       &solve: ' // trim(StringWrite(snes_reason))
    call PrintMsg(option)
    stop_flag = TS_STOP_FAILURE
    return
  endif

  ! print screen output
  call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_FUNCTION, &
                       PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)

  call TimestepperBasePrintStepInfo(this,process_model%output_option, &
                                snes_reason,option)
  write(option%io_buffer,'("  newton = ",i3," [",i8,"]", " linear = ",i5, &
                         &" [",i10,"]")') &
           num_newton_iterations,this%cumulative_newton_iterations, &
           num_linear_iterations,this%cumulative_linear_iterations
  call PrintMsg(option)

  select type(process_model)
    class is (pm_subsurface_flow_type)
      scaled_fnorm = fnorm/process_model%realization_base% &
                     discretization%grid%nmax
    class default
      scaled_fnorm = fnorm
  end select

  option%io_buffer = '  --> SNES Linear/Non-Linear Iterations = ' // &
    trim(StringWrite(num_linear_iterations)) // ' / ' // &
    trim(StringWrite(num_newton_iterations))
  call PrintMsg(option)
  option%io_buffer = '  --> SNES Residual: ' // &
    trim(StringWrite('(e14.6)',fnorm)) // ' ' // &
    trim(StringWrite('(e14.6)',scaled_fnorm)) // ' ' // &
    trim(StringWrite('(e14.6)',inorm))
  call PrintMsg(option)

  option%time = this%target_time
  call process_model%FinalizeTimestep()

  if (this%print_ekg .and. OptionPrintToFile(option)) then
100 format(a32," TIMESTEP ",i10,2es16.8,a,i3,i5,i3,i5,i5,i10)
    write(IUNIT_EKG,100) trim(this%name), this%steps, &
      this%target_time/process_model%output_option%tconv, &
      trim(process_model%output_option%tunit), &
      num_newton_iterations, this%cumulative_newton_iterations, &
      num_linear_iterations, this%cumulative_linear_iterations
  endif

end subroutine TimestepperSteadyStepDT

! ************************************************************************** !

subroutine TimestepperSteadyPrintInfo(this,aux_string,option)
  !
  ! Prints settings for steady timestepper.
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !
  use Timestepper_Base_class
  use Option_module

  implicit none

  class(timestepper_steady_type) :: this
  character(len=*) :: aux_string
  type(option_type) :: option

  call TimestepperSNESPrintInfo(this,'SNES (Steady)',option)
  call SolverPrintNewtonInfo(this%solver,this%name,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)

end subroutine TimestepperSteadyPrintInfo

! ************************************************************************** !

subroutine TimestepperSteadyInputRecord(this)
  !
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/16/20
  !

  implicit none

  class(timestepper_steady_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc steady timestepper: '
  write(id,'(a)') this%name

end subroutine TimestepperSteadyInputRecord

! ************************************************************************** !

recursive subroutine TimestepperSteadyFinalizeRun(this,option)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !
  use Option_module
  use String_module

  implicit none

  class(timestepper_steady_type) :: this
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

#ifdef DEBUG
  call PrintMsg(option,'TimestepperSteadyFinalizeRun()')
#endif

  if (this%cumulative_newton_iterations > 0) then
    string = ' ' // trim(this%name) // &
      ' TS Steady steps = ' // trim(StringWrite(this%steps)) // &
      ' newton = ' // trim(StringWrite(this%cumulative_newton_iterations)) // &
      ' linear = ' // trim(StringWrite(this%cumulative_linear_iterations))
  else if (this%cumulative_linear_iterations > 0) then
    string = ' ' // trim(this%name) // &
      ' TS Steady steps = ' // trim(StringWrite(this%steps)) // &
      ' linear = ' // trim(StringWrite(this%cumulative_linear_iterations))
  else
    string = ' ' // trim(this%name) // &
      ' TS Steady steps = ' // trim(StringWrite(this%steps))
  endif
  call PrintMsg(option,string)
  string = ' ' // trim(this%name) // &
    ' TS Steady time = ' // &
    trim(StringWrite('(f12.1)',this%cumulative_solver_time)) // ' seconds'
  call PrintMsg(option,string)

end subroutine TimestepperSteadyFinalizeRun

! ************************************************************************** !

subroutine TimestepperSteadyDestroy(this)
  !
  ! Deallocates a time stepper
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21
  !

  implicit none

  class(timestepper_steady_type) :: this

  call this%Strip()

!  deallocate(this)
!  nullify(this)

end subroutine TimestepperSteadyDestroy

end module Timestepper_Steady_class
