module Inversion_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module
  use Inversion_TS_Aux_module
  use Option_Inversion_module

  implicit none

  private

  type, public :: inversion_aux_type
    PetscInt :: max_ts
    Mat :: JsensitivityT
    Mat :: M ! solely a pointer
    Vec :: solution ! solely a pointer
    type(inversion_forward_aux_type), pointer :: inversion_forward_aux
    type(inversion_measurement_aux_type), pointer :: measurements(:)
    VecScatter :: scatter_global_to_measurement
    Vec :: measurement_vec
  end type inversion_aux_type

  public :: InversionAuxCreate, &
            InversionAuxDestroy

contains

! ************************************************************************** !

function InversionAuxCreate()
  !
  ! Allocate and initialize auxiliary inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  type(inversion_aux_type), pointer :: InversionAuxCreate

  type(inversion_aux_type), pointer :: aux

  allocate(aux)
  nullify(aux%inversion_forward_aux)

  aux%max_ts = UNINITIALIZED_INTEGER
  aux%M = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC
  nullify(aux%measurements)
  aux%scatter_global_to_measurement = PETSC_NULL_VECSCATTER
  aux%measurement_vec = PETSC_NULL_VEC

  aux%JsensitivityT = PETSC_NULL_MAT

  InversionAuxCreate => aux

end function InversionAuxCreate

! ************************************************************************** !

subroutine InversionAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  type(inversion_aux_type), pointer :: aux

  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  call InversionForwardAuxDestroy(aux%inversion_forward_aux)

  if (aux%JsensitivityT /= PETSC_NULL_MAT) then
    call MatDestroy(aux%JsensitivityT,ierr);CHKERRQ(ierr)
  endif

  ! these objects are destroyed elsewhere, do not destroy
  aux%M = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC
  nullify(aux%measurements)
  aux%scatter_global_to_measurement = PETSC_NULL_VECSCATTER
  aux%measurement_vec = PETSC_NULL_VEC

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
