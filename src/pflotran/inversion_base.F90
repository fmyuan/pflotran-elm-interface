module Inversion_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Driver_module
  use Timer_class

  implicit none

  private

  type, public :: inversion_base_type
    class(driver_type), pointer :: driver
    class(timer_type), pointer :: timer
  contains
    procedure, public :: Init => InversionBaseInit
    procedure, public :: Initialize => InversionBaseInitialize
    procedure, public :: UpdateParameters => InversionBaseUpdateParameters
    procedure, public :: CalculateInverse => InversionBaseCalculateInverse
    procedure, public :: Finalize => InversionBaseFinalize
    procedure, public :: Strip => InversionBaseStrip
  end type inversion_base_type

  public :: InversionBaseInit, &
            InversionBaseFinalize, &
            InversionBaseStrip, &
            InversionBaseDestroy

contains

! ************************************************************************** !

subroutine InversionBaseInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this
  class(driver_type), pointer :: driver

  this%driver => driver
  this%timer => TimerCreate()
  call this%timer%Start()

end subroutine InversionBaseInit

! ************************************************************************** !

subroutine InversionBaseInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseInitialize

! ************************************************************************** !

subroutine InversionBaseUpdateParameters(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseUpdateParameters

! ************************************************************************** !

subroutine InversionBaseCalculateInverse(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseCalculateInverse

! ************************************************************************** !

subroutine InversionBaseFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

  PetscLogDouble :: total_time

  call this%timer%Stop()
  total_time = this%timer%GetCumulativeTime()

end subroutine InversionBaseFinalize

! ************************************************************************** !

subroutine InversionBaseStrip(this)
  !
  ! Deallocates members of inversion base
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

  nullify(this%driver)
  call TimerDestroy(this%timer)

end subroutine InversionBaseStrip

! ************************************************************************** !

subroutine InversionBaseDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionBaseDestroy

end module Inversion_Base_class
