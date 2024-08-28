module Option_Parameter_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: parameter_option_type
    PetscInt :: num_parameters
    character(len=MAXWORDLENGTH), pointer :: parameter_names(:)
  end type parameter_option_type

  public :: OptionParameterCreate, &
            OptionParameterInit, &
            OptionParameterDestroy

contains

! ************************************************************************** !

function OptionParameterCreate()
  !
  ! Allocates and initializes a new OptionParameter object
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/24

  type(parameter_option_type), pointer :: OptionParameterCreate

  type(parameter_option_type), pointer :: option

  allocate(option)

  call OptionParameterInit(option)
  OptionParameterCreate => option

end function OptionParameterCreate

! ************************************************************************** !

subroutine OptionParameterInit(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/24

  type(parameter_option_type) :: option

  option%num_parameters = 0
  nullify(option%parameter_names)

end subroutine OptionParameterInit

! ************************************************************************** !

subroutine OptionParameterDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/24

  type(parameter_option_type), pointer :: option

  if (.not.associated(option)) return

  ! cannot use Utility_module:DeallocateArray due to circular dependency
  if (associated(option%parameter_names)) deallocate(option%parameter_names)
  nullify(option%parameter_names)

  deallocate(option)
  nullify(option)

end subroutine OptionParameterDestroy

end module Option_Parameter_module
