module Factory_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  public :: Initialize, &
            Finalize

contains

! ************************************************************************** !

subroutine Initialize()
  !
  ! Initializes libraries for simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/10/21

  use HDF5_Aux_module
  use Logging_module

  implicit none

  call HDF5Init()
  call LoggingCreate()

end subroutine Initialize

! ************************************************************************** !

subroutine Finalize()
  !
  ! Initializes libraries for simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/10/21
  !
  use HDF5_Aux_module
  use Logging_module

  implicit none

  PetscErrorCode :: ierr

  call LoggingDestroy()
  call HDF5Finalize()
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-options_left','no',ierr);CHKERRQ(ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-objects_left','yes',ierr);CHKERRQ(ierr)
  call PetscFinalize(ierr);CHKERRQ(ierr)

end subroutine Finalize

end module Factory_module
