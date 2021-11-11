module Option_Inversion_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_option_type
  end type inversion_option_type

  public :: OptionInversionCreate, &
            OptionInversionInit, &
            OptionInversionDestroy

contains

! ************************************************************************** !

function OptionInversionCreate()
  !
  ! Allocates and initializes a new OptionInversion object
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/21
  !

  implicit none

  type(inversion_option_type), pointer :: OptionInversionCreate

  type(inversion_option_type), pointer :: option

  allocate(option)

  call OptionInversionInit(option)
  OptionInversionCreate => option

end function OptionInversionCreate

! ************************************************************************** !

subroutine OptionInversionInit(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/21
  !

  implicit none

  type(inversion_option_type) :: option

end subroutine OptionInversionInit

! ************************************************************************** !

subroutine OptionInversionDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 10/14/21
  !

  implicit none

  type(inversion_option_type), pointer :: option

  if (.not.associated(option)) return

  deallocate(option)
  nullify(option)

end subroutine OptionInversionDestroy

end module Option_Inversion_module
