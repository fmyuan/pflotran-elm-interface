module PMC_Geophysics_class

  use Option_module
  use PMC_Base_class
  use PM_ERT_class
  use Realization_Subsurface_class
  use Waypoint_module

  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none


  private

  type, public, extends(pmc_base_type) :: pmc_geophysics_type
    class(realization_subsurface_type), pointer :: realization
    ! current waypoint in pm list (not time_stepper)
    type(waypoint_type), pointer :: cur_waypoint
  contains
    procedure, public :: Init => PMCGeophysicsInit
    procedure, public :: InitializeRun => PMCGeophysicsInitializeRun
    procedure, public :: SetupSolvers => PMCGeophysicsSetupSolvers
    procedure, public :: StepDT => PMCGeophysicsStepDT
    procedure, public :: FinalizeRun => PMCGeophysicsFinalizeRun
    procedure, public :: Destroy => PMCGeophysicsDestroy
  end type pmc_geophysics_type

  public :: PMCGeophysicsCreate, &
            PMCGeophysicsInit

contains

! ************************************************************************** !

function PMCGeophysicsCreate()
  !
  ! Allocates and initializes a new process_model_coupler
  ! object for geophysics.
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
  !

  implicit none

  class(pmc_geophysics_type), pointer :: PMCGeophysicsCreate

  class(pmc_geophysics_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()

  PMCGeophysicsCreate => pmc

end function PMCGeophysicsCreate

! ************************************************************************** !

subroutine PMCGeophysicsInit(this)
  !
  ! Initializes a new process model coupler object.
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
  !
  ! for some reason, Intel with VS want this explicitly specified.
  use PMC_Base_class, only : PMCBaseInit

  implicit none

  class(pmc_geophysics_type) :: this

  call PMCBaseInit(this)
  this%name = 'PMCGeophysics'
  nullify(this%realization)

end subroutine PMCGeophysicsInit

! ************************************************************************** !

recursive subroutine PMCGeophysicsInitializeRun(this)
  !
  ! Last chance to initialize/configure before run
  !
  ! Author: Glenn Hammond
  ! Date: 04/19/21
  !

  implicit none

  class(pmc_geophysics_type) :: this

  class(pm_ert_type), pointer :: pm_ert

  pm_ert => PMERTCast(this%pm_ptr%pm)
  this%cur_waypoint => pm_ert%waypoint_list%first

  ! ensure that the first waypoint is not at time zero.
  if (associated(this%cur_waypoint)) then
    if (this%cur_waypoint%time < 1.d-40) then
      this%option%io_buffer = 'Simulating an ERT survey is not possible at &
        &time zero. Please choose a non-zero time for the first survey.'
      call PrintErrMsg(this%option)
    endif
  endif

  call PMCBaseInitializeRun(this)

end subroutine PMCGeophysicsInitializeRun

! ************************************************************************** !

subroutine PMCGeophysicsSetupSolvers(this)
  !
  ! Author: Glenn Hammond & Piyoosh Jaysaval
  ! Date: 01/29/21
  !
  use Option_module

  implicit none

  class(pmc_geophysics_type) :: this

  type(option_type), pointer :: option
  class(pm_ert_type), pointer :: pm_ert
  PetscErrorCode :: ierr

  option => this%option

  select type(pm=>this%pm_ptr%pm)
    class is(pm_ert_type)
      pm_ert => pm
      ! Sets up solver for ERT
      call pm_ert%SetupSolvers()
    class default
      option%io_buffer = 'Solver setup implemented only for ERT &
                          &geophysics process model.'
      call PrintErrMsg(option)
  end select

end subroutine PMCGeophysicsSetupSolvers

! ************************************************************************** !

subroutine PMCGeophysicsStepDT(this,stop_flag)
  !
  ! Solves a round of steady-state geophysics solutions.
  !
  ! Author: Glenn Hammond & Piyoosh Jaysaval
  ! Date: 01/29/21
  !
  use Option_module
  use Output_Aux_module
  use PM_Base_class
  use Timestepper_Steady_class
  use Timestepper_Base_class, only : TS_STOP_FAILURE
  use Utility_module

  implicit none

  class(pmc_geophysics_type) :: this
  PetscInt :: stop_flag

  class(pm_base_type), pointer :: pm_base
  class(pm_ert_type), pointer :: pm_ert
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(timestepper_steady_type), pointer :: timestepper

  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: local_stop_flag
  PetscInt :: linear_iterations_in_step
  PetscBool :: skip_survey
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE) return

  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)
  call this%PrintHeader()

  option => this%option
  pm_ert => PMERTCast(this%pm_ptr%pm)
  output_option => pm_ert%output_option
  timestepper => TimestepperSteadyCast(this%timestepper)
  linear_iterations_in_step = 0

  skip_survey = PETSC_TRUE
  if (pm_ert%waypoint_list%num_waypoints > 0) then
    if (associated(this%cur_waypoint)) then
      if (Equal(this%cur_waypoint%time,timestepper%target_time)) then
        skip_survey = PETSC_FALSE
        this%cur_waypoint => this%cur_waypoint%next
      endif
    endif
  elseif (option%iflowmode /= NULL_MODE .or. &
          option%itranmode /= NULL_MODE) then
    option%io_buffer = 'SURVEY_TIMES must be listed under &
      &SUBSURFACE_GEOPHYSICS OPTIONS when geophysics is coupled to flow &
      &or transport.'
    call PrintErrMsg(option)
  else
    skip_survey = PETSC_FALSE
  endif

  if (skip_survey) then
    if (this%option%print_screen_flag) then
      write(*, '(/," Time= ",1pe12.5," [",a,"]", &
            &" Skipping geophysics as this is not a survey time.",/)') &
         timestepper%target_time/output_option%tconv,trim(output_option%tunit)
    endif
    if (this%option%print_file_flag) then
      write(this%option%fid_out, '(/," Time= ",1pe12.5," [",a,"]", &
            &" Skipping geophysics as this is not a survey time.",/)') &
         timestepper%target_time/output_option%tconv,trim(output_option%tunit)
    endif
    return
  endif

  call pm_ert%PreSolve()
  call pm_ert%Solve(timestepper%target_time,ierr)
  linear_iterations_in_step = pm_ert%linear_iterations_in_step
  if (ierr /= 0) stop_flag = TS_STOP_FAILURE

  timestepper%steps = timestepper%steps + 1
  timestepper%cumulative_linear_iterations = &
    timestepper%cumulative_linear_iterations + linear_iterations_in_step
  if (this%option%print_screen_flag) then
    write(*, '(/," Step ",i6," Time= ",1pe12.5," [",a,"]", &
         & /,"  linear = ",i5," [",i10,"]")') &
         timestepper%steps, timestepper%target_time/output_option%tconv, &
         trim(output_option%tunit), linear_iterations_in_step, &
         timestepper%cumulative_linear_iterations
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out, '(/," Step ",i6," Time= ",1pe12.5," [",a,"]", &
         & /,"  linear = ",i5," [",i10,"]")') &
         timestepper%steps, this%timestepper%target_time/output_option%tconv, &
         trim(output_option%tunit), linear_iterations_in_step, &
         timestepper%cumulative_linear_iterations
  endif

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
  this%cumulative_time = this%cumulative_time + log_end_time - log_start_time

  ! call this%timer%Start()
  ! do ielectrode = 1, this%num_electrodes
  ! enddo
  ! call this%timer%Stop()

end subroutine PMCGeophysicsStepDT

! ************************************************************************** !

recursive subroutine PMCGeophysicsFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
  !
  use Option_module
  use Timestepper_Steady_class

  implicit none

  class(pmc_geophysics_type) :: this

  class(timestepper_steady_type), pointer :: timestepper
  class(pm_ert_type), pointer :: pm_ert

#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeophysics%FinalizeRun()')
#endif

  nullify(this%realization)
  pm_ert => PMERTCast(this%pm_ptr%pm)
  timestepper => TimestepperSteadyCast(this%timestepper)
  timestepper%cumulative_solver_time = pm_ert%ksp_time
  call PMCBaseFinalizeRun(this)

end subroutine PMCGeophysicsFinalizeRun

! ************************************************************************** !

subroutine PMCGeophysicsStrip(this)
  !
  ! Deallocates members of PMC Geophysics.
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21

  implicit none

  class(pmc_geophysics_type) :: this

  call PMCBaseStrip(this)
  nullify(this%realization)

end subroutine PMCGeophysicsStrip

! ************************************************************************** !

recursive subroutine PMCGeophysicsDestroy(this)
  !
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
  !

  use Option_module

  implicit none

  class(pmc_geophysics_type) :: this

#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeophysics%Destroy()')
#endif

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif

  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif

  call PMCGeophysicsStrip(this)

end subroutine PMCGeophysicsDestroy

! ************************************************************************** !

end module PMC_Geophysics_class
