module Option_Checkpoint_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: CHECKPOINT_BINARY = 1
  PetscInt, parameter, public :: CHECKPOINT_HDF5 = 2
  PetscInt, parameter, public :: CHECKPOINT_BOTH = 3

  type, public :: checkpoint_option_type
    character(len=MAXWORDLENGTH) :: tunit
    PetscReal :: tconv
    PetscReal :: periodic_time_incr
    PetscInt :: periodic_ts_incr
    PetscInt :: format
  end type checkpoint_option_type

  public :: OptionCheckpointCreate, &
            OptionCheckpointInit, &
            OptionCheckpointDestroy

contains

! ************************************************************************** !

function OptionCheckpointCreate()
  !
  ! Allocates and initializes a new OptionCheckpoint object
  !
  ! Author: Glenn Hammond
  ! Date: 11/22/22
  !

  implicit none

  type(checkpoint_option_type), pointer :: OptionCheckpointCreate

  type(checkpoint_option_type), pointer :: option

  allocate(option)

  call OptionCheckpointInit(option)
  OptionCheckpointCreate => option

end function OptionCheckpointCreate

! ************************************************************************** !

subroutine OptionCheckpointInit(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 11/22/22
  !

  implicit none

  type(checkpoint_option_type) :: option

  option%tunit = ''
  option%tconv = 0.d0
  option%periodic_time_incr = UNINITIALIZED_DOUBLE
  option%periodic_ts_incr = 0
  option%format = CHECKPOINT_BINARY

end subroutine OptionCheckpointInit

! ************************************************************************** !

subroutine OptionCheckpointDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 11/22/22
  !

  implicit none

  type(checkpoint_option_type), pointer :: option

  if (.not.associated(option)) return

  deallocate(option)
  nullify(option)

end subroutine OptionCheckpointDestroy

end module Option_Checkpoint_module
