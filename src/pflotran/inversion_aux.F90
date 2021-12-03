module Inversion_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Inversion_TS_Aux_module

  implicit none

  private

  type, public :: inversion_aux_type
    Mat :: JsensitivityT
    Mat :: M ! solely a pointer
    Vec :: solution ! solely a pointer
    PetscInt, pointer :: cell_to_internal_connection(:,:)
    PetscInt, pointer :: cell_to_bc_connection(:,:)
    type(inversion_ts_aux_type), pointer :: inversion_ts_aux_list
  end type inversion_aux_type

  public :: InversionAuxCreate, &
            InversionAuxDestroy

contains

! ************************************************************************** !

function InversionAuxCreate()
  !
  ! Allocate and initialize auxiliary invesion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  type(inversion_aux_type), pointer :: InversionAuxCreate

  type(inversion_aux_type), pointer :: aux

  allocate(aux)
  nullify(aux%cell_to_internal_connection)
  nullify(aux%cell_to_bc_connection)
  nullify(aux%inversion_ts_aux_list)

  aux%M = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC

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

  if (.not.associated(aux)) return

  call DeallocateArray(aux%cell_to_internal_connection)
  call DeallocateArray(aux%cell_to_bc_connection)

  call InversionTSAuxListDestroy(aux%inversion_ts_aux_list,PETSC_TRUE)

  ! these objects are destroyed elsewhere, do not destroy
  aux%JsensitivityT = PETSC_NULL_MAT
  aux%M = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
