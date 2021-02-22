module Timer_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  type, public :: timer_type
    PetscLogDouble :: start_time
    PetscLogDouble :: cumulative_time
  contains
    procedure, public :: Start => TimerStart
    procedure, public :: Stop => TimerStop
    procedure, public :: GetCumulativeTime => TimerGetCumulativeTime
  end type timer_type

  public :: TimerCreate, &
            TimerDestroy

contains

! ************************************************************************** !

function TimerCreate()
  !
  ! Creates a timer object
  !
  ! Author: Glenn Hammond
  ! Date: 01/27/21
  !
  implicit none

  class(timer_type), pointer :: TimerCreate

  class(timer_type), pointer :: timer

  allocate(timer)
  timer%start_time = 0.d0
  timer%cumulative_time = 0.d0

  TimerCreate => timer

end function TimerCreate

! ************************************************************************** !

subroutine TimerStart(this)
  !
  ! Sets start time
  !
  ! Author: Glenn Hammond
  ! Date: 01/27/21
  !
  implicit none

  class(timer_type) :: this

  PetscErrorCode :: ierr

  call PetscTime(this%start_time,ierr);CHKERRQ(ierr)

end subroutine TimerStart

! ************************************************************************** !

subroutine TimerStop(this)
  !
  ! Calculates end time and add difference to cumulative
  !
  ! Author: Glenn Hammond
  ! Date: 01/27/21
  !
  implicit none

  class(timer_type) :: this

  PetscLogDouble :: end_time
  PetscErrorCode :: ierr

  call PetscTime(end_time,ierr);CHKERRQ(ierr)
  this%cumulative_time = this%cumulative_time + (end_time - this%start_time)

end subroutine TimerStop

! ************************************************************************** !

function TimerGetCumulativeTime(this)
  !
  ! Returns cumulative time
  !
  ! Author: Glenn Hammond
  ! Date: 01/27/21
  !
  implicit none

  class(timer_type) :: this

  PetscLogDouble :: TimerGetCumulativeTime

  TimerGetCumulativeTime = this%cumulative_time

end function TimerGetCumulativeTime

! ************************************************************************** !

subroutine TimerDestroy(timer)
  !
  ! Destroys a timer
  !
  ! Author: Glenn Hammond
  ! Date: 01/27/21
  !
  implicit none

  class(timer_type), pointer :: timer

  if (.not.associated(timer)) return

  deallocate(timer)
  nullify(timer)

end subroutine TimerDestroy

end module Timer_class
