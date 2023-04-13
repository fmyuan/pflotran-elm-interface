module Driver_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Communicator_Aux_module
  use PFLOTRAN_Constants_module
  use Print_module

  implicit none

  private

  type, public :: driver_type

    type(print_flags_type), pointer :: print_flags
    type(comm_type), pointer :: comm

    character(len=MAXSTRINGLENGTH) :: input_filename
    character(len=MAXSTRINGLENGTH) :: input_prefix
    character(len=MAXSTRINGLENGTH) :: global_prefix

    PetscInt :: fid_out
    PetscInt :: exit_code                  ! code passed out of PFLOTRAN
                                           ! at end of simulation
    PetscInt :: status

  contains
    procedure, public :: PrintErrMsg => DriverPrintErrorMessage
    procedure, public :: PrintMsg => DriverPrintMessage
    procedure, public :: PrintToScreen => DriverPrintToScreen
    procedure, public :: PrintToFile => DriverPrintToFile
    procedure, public :: IsIORank => DriverIsIORank
  end type driver_type

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
  driver%print_flags => PrintCreateFlags()
  nullify(driver%comm)
  driver%input_filename = ''
  driver%input_prefix = ''
  driver%global_prefix = ''
  driver%fid_out = -1
  driver%exit_code = 0
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
  PetscBool, parameter :: blocking = PETSC_TRUE
  PetscBool, parameter :: byrank = PETSC_FALSE

  call PrintErrorMessage(driver%print_flags,driver%comm,driver%fid_out, &
                         string,driver%exit_code,blocking,byrank)

end subroutine DriverPrintErrorMessage

! ************************************************************************** !

subroutine DriverPrintMessage(driver,string)
  !
  ! Prints a message to the screen and/or file
  !
  ! Author: Glenn Hammond
  ! Date: 11/18/22

  class(driver_type) :: driver
  character(len=*) :: string

  PetscBool, parameter :: advance_ = PETSC_TRUE
  PetscBool, parameter :: byrank = PETSC_FALSE

  call PrintMessage(driver%print_flags,driver%comm,driver%fid_out, &
                    string,advance_,byrank)

end subroutine DriverPrintMessage

! ************************************************************************** !

function DriverPrintToScreen(driver)
  !
  ! Returns boolean indicating whether to print to screen
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type) :: driver

  PetscBool :: DriverPrintToScreen

  DriverPrintToScreen = PrintToScreen(driver%print_flags,driver%comm)

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

  DriverPrintToFile = PrintToFile(driver%print_flags,driver%comm)

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

  DriverIsIORank = CommIsIORank(driver%comm)

end function DriverIsIORank

! ************************************************************************** !

function DriverCreatePrintHandler(driver)
  !
  ! Creates a print handler with flags for file and screen io
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23
  !
  use Print_module

  implicit none

  class(driver_type) :: driver

  type(print_handler_type), pointer :: DriverCreatePrintHandler

  DriverCreatePrintHandler => &
    PrintCreateHandler(driver%print_flags, &
                       driver%comm, &
                       driver%fid_out, &
                       driver%exit_code, &
                       PETSC_TRUE, & ! advance
                       PETSC_TRUE, & ! blocking
                       PETSC_FALSE)  ! byrank

end function DriverCreatePrintHandler

! ************************************************************************** !

subroutine DriverDestroy(driver)
  !
  ! Deallocates an driver
  !
  ! Author: Glenn Hammond
  ! Date: 05/26/21

  class(driver_type), pointer :: driver

  call PrintDestroyFlags(driver%print_flags)
  call CommDestroy(driver%comm)
  if (driver%fid_out > 0) close(driver%fid_out)

  deallocate(driver)
  nullify(driver)

end subroutine DriverDestroy

end module Driver_class
