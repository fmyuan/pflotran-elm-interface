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
    procedure, public :: CheckpointBinary => PMCGeophysicsCheckpointBinary
    procedure, public :: RestartBinary => PMCGeophysicsRestartBinary
    procedure, public :: CheckpointHDF5 => PMCGeophysicsCheckpointHDF5
    procedure, public :: RestartHDF5 => PMCGeophysicsRestartHDF5
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

  ! if restarting at a time greater than time 0, the waypoint pointer needs
  ! to skip ahead
  if (this%option%restart_flag .and. &
      .not. pm_ert%skip_restart .and. &
      .not.Initialized(this%option%restart_time)) then
    call WaypointSkipToTime(this%cur_waypoint, &
                            this%timestepper%target_time)
  endif

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
  use Option_Inversion_module
  use Output_Aux_module
  use PM_Base_class
  use Timestepper_Steady_class
  use Timestepper_Base_class, only : TS_STOP_FAILURE
  use Utility_module

  implicit none

  class(pmc_geophysics_type) :: this
  PetscInt :: stop_flag

  class(pm_ert_type), pointer :: pm_ert
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(timestepper_steady_type), pointer :: timestepper

  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
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

  if (associated(option%inversion)) then
    if (option%inversion%coupled_flow_ert .and. &
        .not.option%inversion%calculate_ert) then
      write(option%io_buffer,'(" Time= ",1pe12.5," [",a,"]", &
            &" Skipping geophysics as this is a perturbed coupled &
            &flow-ert run.")') &
        timestepper%target_time/output_option%tconv,trim(output_option%tunit)
      call PrintMsg(option)
      return
    endif
  endif

  if (skip_survey) then
    write(this%option%io_buffer,'(" Time= ",1pe12.5," [",a,"]", &
          &" Skipping geophysics as this is not a survey time.")') &
      timestepper%target_time/output_option%tconv,trim(output_option%tunit)
    call PrintMsg(this%option)
    return
  endif

  call pm_ert%PreSolve()
  call pm_ert%Solve(timestepper%target_time,ierr)
  linear_iterations_in_step = pm_ert%linear_iterations_in_step
  if (ierr /= 0) stop_flag = TS_STOP_FAILURE

  timestepper%steps = timestepper%steps + 1
  timestepper%cumulative_linear_iterations = &
    timestepper%cumulative_linear_iterations + linear_iterations_in_step
  write(this%option%io_buffer,'(" Step ",i6," Time= ",1pe12.5," [",a,"]", &
                              &a,"  linear = ",i5," [",i10,"]")') &
       timestepper%steps, timestepper%target_time/output_option%tconv, &
       trim(output_option%tunit),new_line('a'), &
       linear_iterations_in_step,timestepper%cumulative_linear_iterations
  call PrintMsg(this%option)

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

recursive subroutine PMCGeophysicsCheckpointBinary(this,viewer,append_name)
  !
  ! Checkpoints geophysics PMC, timestepper and state variables.
  !
  ! Author: Glenn Hammond
  ! Date: 09/13/22
  !
  use PM_Base_class

  implicit none

  class(pmc_geophysics_type) :: this
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: append_name

  class(pm_base_type), pointer :: cur_pm

  ! if the top PMC
  if (this%is_master) then
    this%option%io_buffer = 'PMC Geophysics cannot checkpoint as the master &
      &process model coupler'
    call PrintErrMsg(this%option)
  endif

  if (associated(this%timestepper)) then
    call this%timestepper%CheckpointBinary(viewer,this%option)
  endif

  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%CheckpointBinary(viewer)
    cur_pm => cur_pm%next
  enddo

  if (associated(this%child)) then
    call this%child%CheckpointBinary(viewer,append_name)
  endif

  if (associated(this%peer)) then
    call this%peer%CheckpointBinary(viewer,append_name)
  endif

end subroutine PMCGeophysicsCheckpointBinary

! ************************************************************************** !

recursive subroutine PMCGeophysicsRestartBinary(this,viewer)
  !
  ! Restarts PMC timestepper and state variables.
  !
  ! Author: Glenn Hammond
  ! Date: 09/13/22
  !
  use PM_Base_class

  implicit none

  class(pmc_geophysics_type) :: this
  PetscViewer :: viewer

  class(pm_base_type), pointer :: cur_pm

  ! if the top PMC
  if (this%is_master) then
    this%option%io_buffer = 'PMC Geophysics cannot restart as the master &
      &process model coupler'
    call PrintErrMsg(this%option)
  endif

  if (associated(this%timestepper)) then
    call this%timestepper%RestartBinary(viewer,this%option)
    if (Initialized(this%option%restart_time)) then
      ! simply a flag to set time back to zero, no matter what the restart
      ! time is set to.
      call this%timestepper%Reset()
      ! note that this sets the target time back to zero.
    endif
! this skip occurs later
!    call WaypointSkipToTime(this%cur_waypoint, &
!                            this%timestepper%target_time)
    this%option%time = this%timestepper%target_time
  endif

  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    if (cur_pm%skip_restart) then
      this%option%io_buffer = 'Due to sequential nature of binary files, &
        &skipping restart for binary formatted files is not allowed.'
      call PrintErrMsg(this%option)
    endif
    call cur_pm%RestartBinary(viewer)
    cur_pm => cur_pm%next
  enddo

  if (associated(this%child)) then
    call this%child%RestartBinary(viewer)
  endif

  if (associated(this%peer)) then
    call this%peer%RestartBinary(viewer)
  endif

end subroutine PMCGeophysicsRestartBinary

! ************************************************************************** !

recursive subroutine PMCGeophysicsCheckpointHDF5(this,h5_chk_grp_id,append_name)
  !
  ! Checkpoints PMC timestepper and state variables in HDF5 format.
  !
  ! Author: Glenn Hammond
  ! Date: 09/13/22
  !
  use hdf5
  use HDF5_Aux_module
  use PM_Base_class

  implicit none

  class(pmc_geophysics_type) :: this
  integer(HID_T) :: h5_chk_grp_id
  character(len=MAXSTRINGLENGTH) :: append_name

  integer(HID_T) :: h5_pmc_grp_id
  integer(HID_T) :: h5_pm_grp_id

  class(pm_base_type), pointer :: cur_pm

  ! if the top PMC
  if (this%is_master) then
    this%option%io_buffer = 'PMC Geophysics cannot checkpoint as the master &
      &process model coupler'
    call PrintErrMsg(this%option)
  else
    call HDF5GroupCreate(h5_chk_grp_id, trim(this%name), &
                         h5_pmc_grp_id, this%option)
  endif

  if (associated(this%timestepper)) then
    call this%timestepper%CheckpointHDF5(h5_pmc_grp_id, this%option)
  endif

  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit

    call HDF5GroupCreate(h5_pmc_grp_id, trim(cur_pm%name), h5_pm_grp_id, &
                         this%option)
    call cur_pm%CheckpointHDF5(h5_pm_grp_id)
    call HDF5GroupClose(h5_pm_grp_id,this%option)

    cur_pm => cur_pm%next
  enddo

  call HDF5GroupClose(h5_pmc_grp_id,this%option)

  if (associated(this%child)) then
    call this%child%CheckpointHDF5(h5_chk_grp_id,append_name)
  endif

  if (associated(this%peer)) then
    call this%peer%CheckpointHDF5(h5_chk_grp_id,append_name)
  endif

end subroutine PMCGeophysicsCheckpointHDF5

! ************************************************************************** !

recursive subroutine PMCGeophysicsRestartHDF5(this,h5_chk_grp_id)
  !
  ! Restarts PMC timestepper and state variables from a HDF5
  !
  ! Author: Glenn Hammond
  ! Date: 09/13/22
  !
  use hdf5
  use HDF5_Aux_module
  use PM_Base_class

  implicit none

  class(pmc_geophysics_type) :: this
  integer(HID_T) :: h5_chk_grp_id

  class(pm_base_type), pointer :: cur_pm

  integer(HID_T) :: h5_pmc_grp_id
  integer(HID_T) :: h5_pm_grp_id

  PetscBool :: skip_restart

  ! search pm for skip restart flag which will apply to everything in pmc
  skip_restart = PETSC_FALSE
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    if (cur_pm%skip_restart) then
      skip_restart = PETSC_TRUE
      exit
    endif
    cur_pm => cur_pm%next
  enddo

  ! if the top PMC
  if (this%is_master) then
    this%option%io_buffer = 'PMC Geophysics cannot restart as the master &
      &process model coupler'
    call PrintErrMsg(this%option)
  else
    if (.not.skip_restart) then
      call HDF5GroupOpen(h5_chk_grp_id,this%name,h5_pmc_grp_id, &
                         this%option%driver)
    endif
  endif

  if (associated(this%timestepper)) then
    if (.not.skip_restart) then
      call this%timestepper%RestartHDF5(h5_pmc_grp_id, this%option)
    endif

    if (Initialized(this%option%restart_time)) then
      ! simply a flag to set time back to zero, no matter what the restart
      ! time is set to.
      call this%timestepper%Reset()
      ! note that this sets the target time back to zero.
    else if (skip_restart) then
        this%option%io_buffer = 'Restarted simulations that SKIP_RESTART on &
          &checkpointed process models must restart at time 0.'
        call PrintErrMsg(this%option)
    endif
! this skip occurs later
!    call WaypointSkipToTime(this%cur_waypoint, &
!                            this%timestepper%target_time)
    this%option%time = this%timestepper%target_time
  endif

  if (.not.skip_restart) then
    cur_pm => this%pm_list
    do
      if (.not.associated(cur_pm)) exit
      call HDF5GroupOpen(h5_pmc_grp_id,cur_pm%name,h5_pm_grp_id, &
                         this%option%driver)
      call cur_pm%RestartHDF5(h5_pm_grp_id)
      call HDF5GroupClose(h5_pm_grp_id,this%option)
      cur_pm => cur_pm%next
    enddo
    call HDF5GroupClose(h5_pmc_grp_id,this%option)
  endif

  if (associated(this%child)) then
    call this%child%RestartHDF5(h5_chk_grp_id)
  endif

  if (associated(this%peer)) then
    call this%peer%RestartHDF5(h5_chk_grp_id)
  endif

end subroutine PMCGeophysicsRestartHDF5

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
