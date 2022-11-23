module Driver_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Communicator_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: driver_type

    type(comm_type), pointer :: comm

    character(len=MAXSTRINGLENGTH) :: input_filename
    character(len=MAXSTRINGLENGTH) :: input_prefix
    character(len=MAXSTRINGLENGTH) :: global_prefix

    PetscInt :: fid_out
    PetscBool :: print_to_screen
    PetscBool :: print_to_file
    PetscMPIInt :: io_rank
    PetscInt :: exit_code                  ! code passed out of PFLOTRAN
                                           ! at end of simulation
    PetscInt :: status

  contains
    procedure, public :: PrintErrMsg => DriverPrintErrorMessage
    procedure, public :: PrintMsg => DriverPrintMessage1
    procedure, public :: PrintToScreen => DriverPrintToScreen
    procedure, public :: PrintToFile => DriverPrintToFile
    procedure, public :: IsIORank => DriverIsIORank
  end type driver_type

  interface DriverPrintMessage
    module procedure :: DriverPrintMessage1
    module procedure :: DriverPrintMessage2
    module procedure :: DriverPrintMessage3
  end interface

  public :: DriverPrintMessage

  public :: DriverCreate, &
            DriverSetComm, &
            DriverDestroy

contains

! ************************************************************************** !

function DriverCreate()
  !
  ! Allocates and initializes a new Driver object
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type), pointer :: DriverCreate

  class(driver_type), pointer :: driver

  allocate(driver)
  nullify(driver%comm)
  driver%input_filename = ''
  driver%input_prefix = ''
  driver%global_prefix = ''
  driver%fid_out = -1
  driver%print_to_screen = PETSC_TRUE
  driver%print_to_file = PETSC_TRUE
  driver%exit_code = 0
  driver%io_rank = 0
  driver%status = 0

  DriverCreate => driver

end function DriverCreate

! ************************************************************************** !

subroutine DriverSetComm(driver,comm)
  !
  ! Sets comm pointer
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type) :: driver
  type(comm_type), pointer :: comm

  driver%comm => comm

end subroutine DriverSetComm

! ************************************************************************** !

subroutine DriverPrintErrorMessage(driver,string)
  !
  ! Prints an error message
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type) :: driver
  character(len=*) :: string

  PetscBool :: petsc_initialized
  PetscErrorCode :: ierr

  if (driver%PrintToScreen()) then
    print *
    print *, 'ERROR: ' // trim(string)
    print *
    print *, 'Stopping!'
  endif
  call MPI_Barrier(driver%comm%mycomm,ierr);CHKERRQ(ierr)
  call PetscInitialized(petsc_initialized,ierr);CHKERRQ(ierr)
  if (petsc_initialized) then
    call PetscFinalize(ierr);CHKERRQ(ierr)
  endif
  select case(driver%exit_code)
    case(EXIT_FAILURE)
      call exit(driver%exit_code)
    case default
      call exit(EXIT_USER_ERROR)
  end select

end subroutine DriverPrintErrorMessage

! ************************************************************************** !

subroutine DriverPrintMessage1(driver,string)
  !
  ! Prints a message to the screen and/or file
  !
  ! Author: Glenn Hammond
  ! Date: 11/18/22

  class(driver_type) :: driver
  character(len=*) :: string

  PetscBool, parameter :: advance_ = PETSC_TRUE

  call DriverPrintMessage3(driver,driver%fid_out,string,advance_)

end subroutine DriverPrintMessage1

! ************************************************************************** !

subroutine DriverPrintMessage2(driver,string,advance_)
  !
  ! Prints a message to the screen and/or file
  !
  ! Author: Glenn Hammond
  ! Date: 11/18/22

  class(driver_type) :: driver
  character(len=*) :: string
  PetscBool :: advance_

  call DriverPrintMessage3(driver,driver%fid_out,string,advance_)

end subroutine DriverPrintMessage2

! ************************************************************************** !

subroutine DriverPrintMessage3(driver,fid,string,advance_)
  !
  ! Prints a message to the screen and/or file
  !
  ! Author: Glenn Hammond
  ! Date: 11/18/22

  class(driver_type) :: driver
  PetscInt :: fid
  character(len=*) :: string
  PetscBool :: advance_

  if (driver%PrintToScreen()) then
    if (.not.advance_) then
      write(STDOUT_UNIT,'(a)',advance='no') trim(string)
    else
      write(STDOUT_UNIT,'(a)') trim(string)
    endif
  endif
  if (driver%PrintToFile() .and. fid > 0) then
    if (.not.advance_) then
      write(fid,'(a)',advance='no') trim(string)
    else
      write(fid,'(a)') trim(string)
    endif
  endif

end subroutine DriverPrintMessage3

! ************************************************************************** !

function DriverPrintToScreen(driver)
  !
  ! Returns boolean indicating whether to print to screen
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type) :: driver

  PetscBool :: DriverPrintToScreen

  if (driver%comm%myrank == driver%io_rank .and. driver%print_to_screen) then
    DriverPrintToScreen = PETSC_TRUE
  else
    DriverPrintToScreen = PETSC_FALSE
  endif

end function DriverPrintToScreen

! ************************************************************************** !

function DriverPrintToFile(driver)
  !
  ! Returns boolean indicating whether to print to file
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type) :: driver

  PetscBool :: DriverPrintToFile

  if (driver%comm%myrank == driver%io_rank .and. driver%print_to_file) then
    DriverPrintToFile = PETSC_TRUE
  else
    DriverPrintToFile = PETSC_FALSE
  endif

end function DriverPrintToFile

! ************************************************************************** !

function DriverIsIORank(driver)
  !
  ! Returns boolean indicating whether to print to file
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type) :: driver

  PetscBool :: DriverIsIORank

  DriverIsIORank = (driver%comm%myrank == driver%io_rank)

end function DriverIsIORank

! ************************************************************************** !

subroutine DriverDestroy(driver)
  !
  ! Deallocates an driver
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type), pointer :: driver

  call CommDestroy(driver%comm)
  if (driver%fid_out > 0) close(driver%fid_out)

  deallocate(driver)
  nullify(driver)

end subroutine DriverDestroy

end module Driver_class
