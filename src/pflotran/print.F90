module Print_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Communicator_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: print_flags_type
    PetscBool :: print_to_screen
    PetscBool :: print_to_file
  end type print_flags_type

  interface PrintCreateFlags
    module procedure :: PrintCreateFlags1
    module procedure :: PrintCreateFlags2
  end interface PrintCreateFlags

  interface PrintInitFlags
    module procedure :: PrintInitFlags1
    module procedure :: PrintInitFlags2
  end interface PrintInitFlags

  public :: PrintCreateFlags, &
            PrintInitFlags, &
            PrintToScreen, &
            PrintToFile, &
            PrintSetPrintToScreenFlag, &
            PrintSetPrintToFileFlag, &
            PrintMessage, &
            PrintErrorMessage, &
            PrintDestroyFlags


contains

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

subroutine PrintErrorMessage(print_flags,comm,fid,string,exit_code, &
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
    local_string = 'ERROR(' // StringWrite(comm%rank) // '): ' // trim(string)
    if (print_flags%print_to_screen) then
      write(STDOUT_UNIT,'(/,a,//,"Stopping!")') local_string
    endif
    if (print_flags%print_to_file .and. fid > 0) then
      write(fid,'(/,a,//,"Stopping!")') local_string
    endif
  else
    local_string = 'ERROR: ' // trim(string)
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

end subroutine PrintErrorMessage

! ************************************************************************** !

subroutine PrintMessage(print_flags,comm,fid,string,advance_,byrank)
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
    local_string = '(' // StringWrite(comm%rank) // '): ' // trim(string)
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

end subroutine PrintMessage

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
