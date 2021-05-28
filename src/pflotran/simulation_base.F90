module Simulation_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Driver_module
  use Option_module
  use Output_Aux_module
  use Output_module
  use Simulation_Aux_module
  use Timer_class
  use Waypoint_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: simulation_base_type
    class(driver_type), pointer :: driver
    type(option_type), pointer :: option
    type(waypoint_list_type), pointer :: waypoint_list_outer ! for outer sync loop
    class(timer_type), pointer :: timer
  contains
    procedure, public :: Init => SimulationBaseInit
    procedure, public :: InitializeRun => SimulationBaseInitializeRun
    procedure, public :: InputRecord => SimulationBaseInputRecord
    procedure, public :: ExecuteRun => SimulationBaseExecuteRun
    procedure, public :: RunToTime => SimulationBaseRunToTime
    procedure, public :: FinalizeRun => SimulationBaseFinalizeRun
    procedure, public :: Strip => SimulationBaseStrip
  end type simulation_base_type

  public :: SimulationBaseCreate, &
            SimulationBaseInit, &
            SimulationBaseInitializeRun, &
            SimulationBaseInputRecord, &
            SimulationBaseFinalizeRun, &
            SimulationBaseStrip, &
            SimulationBaseDestroy

  public :: SimulationBaseInputRecordPrint

contains

! ************************************************************************** !

function SimulationBaseCreate(driver,option)
  !
  ! Allocates and initializes a new simulation object
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !

  use Option_module

  implicit none

  class(simulation_base_type), pointer :: SimulationBaseCreate

  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  allocate(SimulationBaseCreate)
  call SimulationBaseCreate%Init(driver,option)

end function SimulationBaseCreate

! ************************************************************************** !

subroutine SimulationBaseInit(this,driver,option)
  !
  ! Initializes a new simulation object
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Option_module
  use Waypoint_module

  implicit none

  class(simulation_base_type) :: this
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  this%driver => driver
  this%option => option
  this%waypoint_list_outer => WaypointListCreate()
  this%timer => TimerCreate()

end subroutine SimulationBaseInit

! ************************************************************************** !

subroutine SimulationBaseInitializeRun(this)
  !
  ! Initializes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Option_module

  implicit none

  class(simulation_base_type) :: this

  call this%timer%Start()

end subroutine SimulationBaseInitializeRun

! ************************************************************************** !

subroutine SimulationBaseInputRecordPrint(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !
  implicit none

  class(simulation_base_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  PetscBool :: is_open

  inquire(id, OPENED=is_open)
  if (is_open .and. OptionPrintToFile(this%option)) then
  !----------------------------------------------------------------------------
    if (OptionPrintToScreen(this%option)) then
      write (*,*) 'Printing input record file.'
    endif

    ! print simulation-specific information
    call this%InputRecord()
  !----------------------------------------------------------------------------
  endif

end subroutine SimulationBaseInputRecordPrint

! ************************************************************************** !

subroutine SimulationBaseInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  ! This subroutine must be extended in the extended simulation objects.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !
  implicit none

  class(simulation_base_type) :: this

  this%option%io_buffer = 'SimulationBaseInputRecord must be extended for &
    &each simulation mode.'
  call PrintErrMsg(this%option)

end subroutine SimulationBaseInputRecord

! ************************************************************************** !

subroutine SimulationBaseExecuteRun(this)
  !
  ! Initializes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  implicit none

  class(simulation_base_type) :: this

  this%option%io_buffer = 'SimulationExecuteRun must be extended for &
    &each simulation mode.'
  call PrintErrMsg(this%option)

end subroutine SimulationBaseExecuteRun

! ************************************************************************** !

subroutine SimulationBaseRunToTime(this,target_time)
  !
  ! Executes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Option_module

  implicit none

  class(simulation_base_type) :: this
  PetscReal :: target_time

  this%option%io_buffer = 'SimulationExecuteRun must be extended for &
    &each simulation mode.'
  call PrintErrMsg(this%option)

end subroutine SimulationBaseRunToTime

! ************************************************************************** !

subroutine SimulationBaseFinalizeRun(this)
  !
  ! Finalizes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !


  implicit none

  class(simulation_base_type) :: this

  PetscLogDouble :: total_time

  call this%timer%Stop()
  total_time = this%timer%GetCumulativeTime()

  if (this%option%myrank == this%option%io_rank) then

    if (this%option%print_to_screen) then
      write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        total_time, &
        total_time/60.d0, &
        total_time/3600.d0
    endif
    if (this%option%print_to_file) then
      write(this%option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        total_time, &
        total_time/60.d0, &
        total_time/3600.d0
    endif
  endif

end subroutine SimulationBaseFinalizeRun

! ************************************************************************** !

subroutine SimulationBaseStrip(this)
  !
  ! Deallocates members of simulation base
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  implicit none

  class(simulation_base_type) :: this

  nullify(this%driver)
  nullify(this%option)
  call WaypointListDestroy(this%waypoint_list_outer)
  call TimerDestroy(this%timer)

end subroutine SimulationBaseStrip

! ************************************************************************** !

subroutine SimulationBaseDestroy(simulation)
  !
  ! Deallocates a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  implicit none

  class(simulation_base_type), pointer :: simulation

  if (.not.associated(simulation)) return

  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)

end subroutine SimulationBaseDestroy

end module Simulation_Base_class
