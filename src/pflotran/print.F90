module Print_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Communicator_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: print_handler_type
    type(print_flags_type), pointer :: print_flags
    type(comm_type), pointer :: comm
    PetscInt :: fid
    PetscInt :: exit_code
    PetscBool :: advance
    PetscBool :: blocking
    PetscBool :: byrank
  end type print_handler_type

  type, public :: print_flags_type
    PetscBool :: print_to_screen
    PetscBool :: print_to_file
  end type print_flags_type

  interface PrintCreateHandler
    module procedure :: PrintCreateHandler1
    module procedure :: PrintCreateHandler2
  end interface PrintCreateHandler


  interface PrintCreateFlags
    module procedure :: PrintCreateFlags1
    module procedure :: PrintCreateFlags2
  end interface PrintCreateFlags

  interface PrintInitFlags
    module procedure :: PrintInitFlags1
    module procedure :: PrintInitFlags2
  end interface PrintInitFlags

  interface PrintErrorMessage
    module procedure :: PrintErrorMessage1
    module procedure :: PrintErrorMessage2
  end interface PrintErrorMessage

  interface PrintMessage
    module procedure :: PrintMessage1
    module procedure :: PrintMessage2
  end interface PrintMessage

  public :: PrintCreateHandler, &
            PrintCreateFlags, &
            PrintInitFlags, &
            PrintToScreen, &
            PrintToFile, &
            PrintSetPrintToScreenFlag, &
            PrintSetPrintToFileFlag, &
            PrintMessage, &
            PrintErrorMessage, &
            PrintDestroyFlags, &
            PrintDestroyHandler

contains

! ************************************************************************** !

function PrintCreateHandler1()
  !
  ! Creates the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23

  type(print_handler_type), pointer :: PrintCreateHandler1

  type(print_handler_type), pointer :: handler

  allocate(handler)
  handler%fid = UNINITIALIZED_INTEGER
  handler%exit_code = UNINITIALIZED_INTEGER
  handler%advance = PETSC_FALSE
  handler%blocking = PETSC_FALSE
  handler%byrank = PETSC_FALSE
  nullify(handler%print_flags)
  nullify(handler%comm)

  PrintCreateHandler1 => handler

end function PrintCreateHandler1

! ************************************************************************** !

function PrintCreateHandler2(print_flags,comm,fid,exit_code,advance, &
                             blocking,byrank)
  !
  ! Creates the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23

  type(print_flags_type), pointer :: print_flags
  type(comm_type), pointer :: comm
  PetscInt :: fid
  PetscInt :: exit_code
  PetscBool :: advance
  PetscBool :: blocking
  PetscBool :: byrank

  type(print_handler_type), pointer :: PrintCreateHandler2

  type(print_handler_type), pointer :: handler

  handler => PrintCreateHandler()
  handler%print_flags => print_flags
  handler%comm => comm
  handler%fid = fid
  handler%exit_code = exit_code
  handler%advance = advance
  handler%blocking = blocking
  handler%byrank = byrank

  PrintCreateHandler2 => handler

end function PrintCreateHandler2

! ************************************************************************** !

function PrintCreateFlags1()
  !
  ! Creates the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type), pointer :: PrintCreateFlags1

  allocate(PrintCreateFlags1)
  call PrintInitFlags(PrintCreateFlags1)

end function PrintCreateFlags1

! ************************************************************************** !

function PrintCreateFlags2(print_flags)
  !
  ! Creates the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags

  type(print_flags_type), pointer :: PrintCreateFlags2

  allocate(PrintCreateFlags2)
  call PrintInitFlags(PrintCreateFlags2,print_flags)

end function PrintCreateFlags2

! ************************************************************************** !

subroutine PrintInitFlags1(print_flags)
  !
  ! Initializes the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags

  print_flags%print_to_screen = PETSC_TRUE
  print_flags%print_to_file = PETSC_TRUE

end subroutine PrintInitFlags1

! ************************************************************************** !

subroutine PrintInitFlags2(print_flags,older_print_flags)
  !
  ! Initializes the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags
  type(print_flags_type) :: older_print_flags

  print_flags%print_to_screen = older_print_flags%print_to_screen
  print_flags%print_to_file = older_print_flags%print_to_file

end subroutine PrintInitFlags2

! ************************************************************************** !

function PrintToScreen(print_flags,comm)
  !
  ! Toggles printing to the screen
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags
  type(comm_type) :: comm

  PetscBool :: PrintToScreen

  PrintToScreen = CommIsIORank(comm) .and. print_flags%print_to_screen

end function PrintToScreen

! ************************************************************************** !

function PrintToFile(print_flags,comm)
  !
  ! Toggles printing to a file
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags
  type(comm_type) :: comm

  PetscBool :: PrintToFile

  PrintToFile = CommIsIORank(comm) .and. print_flags%print_to_file

end function PrintToFile

! ************************************************************************** !

subroutine PrintSetPrintToScreenFlag(print_flags,flag)
  !
  ! Sets the toggle for printing to the screen
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags
  PetscBool :: flag

  print_flags%print_to_screen = flag

end subroutine PrintSetPrintToScreenFlag

! ************************************************************************** !

subroutine PrintSetPrintToFileFlag(print_flags,flag)
  !
  ! Sets the toggle for printing to the file
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type) :: print_flags
  PetscBool :: flag

  print_flags%print_to_file = flag

end subroutine PrintSetPrintToFileFlag

! ************************************************************************** !

subroutine PrintErrorMessage1(print_handler,string)
  !
  ! Prints an error message
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23

  type(print_handler_type) :: print_handler
  character(len=*) :: string

  call PrintErrorMessage(print_handler%print_flags, &
                         print_handler%comm, &
                         print_handler%fid, &
                         string, &
                         print_handler%exit_code, &
                         print_handler%blocking, &
                         print_handler%byrank)

end subroutine PrintErrorMessage1

! ************************************************************************** !

subroutine PrintErrorMessage2(print_flags,comm,fid,string,exit_code, &
                              blocking,byrank)
  !
  ! Prints an error message
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  use String_module

  type(print_flags_type) :: print_flags
  type(comm_type) :: comm
  PetscInt :: fid
  character(len=*) :: string
  PetscInt :: exit_code
  PetscBool :: blocking
  PetscBool :: byrank

  PetscBool :: petsc_initialized
  character(len=:), allocatable :: local_string
  PetscErrorCode :: ierr

  if (byrank) then
    allocate(local_string,source = 'ERROR(' //  StringWrite(comm%rank) // &
                                   '): ' // trim(string))
    if (print_flags%print_to_screen) then
      write(STDOUT_UNIT,'(/,a,//,"Stopping!")') local_string
    endif
    if (print_flags%print_to_file .and. fid > 0) then
      write(fid,'(/,a,//,"Stopping!")') local_string
    endif
  else
    allocate(local_string,source = 'ERROR: ' // trim(string))
    if (PrintToScreen(print_flags,comm)) then
      write(STDOUT_UNIT,'(/,a,//,"Stopping!")') local_string
    endif
    if (PrintToFile(print_flags,comm) .and. fid > 0) then
      write(fid,'(/,a,//,"Stopping!")') local_string
    endif
  endif

  if (blocking) then
    call MPI_Barrier(comm%communicator,ierr);CHKERRQ(ierr)
    call PetscInitialized(petsc_initialized,ierr);CHKERRQ(ierr)
    if (petsc_initialized) then
      call PetscFinalize(ierr);CHKERRQ(ierr)
    endif
    select case(exit_code)
      case(EXIT_FAILURE)
        call exit(exit_code)
      case default
        call exit(EXIT_USER_ERROR)
    end select
  endif

end subroutine PrintErrorMessage2

! ************************************************************************** !

subroutine PrintMessage1(print_handler,string)
  !
  ! Prints a message
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23

  type(print_handler_type) :: print_handler
  character(len=*) :: string

  call PrintMessage(print_handler%print_flags, &
                    print_handler%comm, &
                    print_handler%fid, &
                    string, &
                    print_handler%advance, &
                    print_handler%byrank)

end subroutine PrintMessage1

! ************************************************************************** !

subroutine PrintMessage2(print_flags,comm,fid,string,advance_,byrank)
  !
  ! Prints a message to the screen and/or file
  !
  ! Author: Glenn Hammond
  ! Date: 02/17/23

  use String_module

  type(print_flags_type) :: print_flags
  type(comm_type) :: comm
  PetscInt :: fid
  character(len=*) :: string
  PetscBool :: advance_
  PetscBool :: byrank

  character(len=:), allocatable :: local_string

  if (byrank) then
    allocate(local_string,source = '(' // StringWrite(comm%rank) // &
                                   '): ' // trim(string))
    if (print_flags%print_to_screen) then
      if (.not.advance_) then
        write(STDOUT_UNIT,'(a)',advance='no') trim(local_string)
      else
        write(STDOUT_UNIT,'(a)') trim(local_string)
      endif
    endif
    if (print_flags%print_to_file .and. fid > 0) then
      if (.not.advance_) then
        write(fid,'(a)',advance='no') trim(local_string)
      else
        write(fid,'(a)') trim(local_string)
      endif
    endif
  else
    if (PrintToScreen(print_flags,comm)) then
      if (.not.advance_) then
        write(STDOUT_UNIT,'(a)',advance='no') trim(string)
      else
        write(STDOUT_UNIT,'(a)') trim(string)
      endif
    endif
    if (PrintToFile(print_flags,comm) .and. fid > 0) then
      if (.not.advance_) then
        write(fid,'(a)',advance='no') trim(string)
      else
        write(fid,'(a)') trim(string)
      endif
    endif
  endif

end subroutine PrintMessage2

! ************************************************************************** !

subroutine PrintDestroyHandler(print_handler)
  !
  ! Destroys the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 04/13/23

  type(print_handler_type), pointer :: print_handler

  if (.not.associated(print_handler)) return

  deallocate(print_handler)
  nullify(print_handler)

end subroutine PrintDestroyHandler

! ************************************************************************** !

subroutine PrintDestroyFlags(print_flags)
  !
  ! Destroys the print flags object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/23

  type(print_flags_type), pointer :: print_flags

  if (.not.associated(print_flags)) return

  deallocate(print_flags)
  nullify(print_flags)

end subroutine PrintDestroyFlags

end module Print_module
