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

  type(inversion_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_ts_aux_type), pointer :: next_inversion_ts_aux

  if (.not.associated(aux)) return

  call DeallocateArray(aux%cell_to_internal_connection)
  call DeallocateArray(aux%cell_to_bc_connection)

  cur_inversion_ts_aux => aux%inversion_ts_aux_list
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    next_inversion_ts_aux => cur_inversion_ts_aux%next
    if (.not.associated(next_inversion_ts_aux)) then
      print *
      print *, '  Last Inversion TS Aux timestep and time: ', &
        cur_inversion_ts_aux%timestep, cur_inversion_ts_aux%time
      print *
    endif
    call InversionTSAuxDestroy(cur_inversion_ts_aux)
    cur_inversion_ts_aux => next_inversion_ts_aux
  enddo

  ! these objects are destroyed elsewhere, do not destroy
  aux%JsensitivityT = PETSC_NULL_MAT
  aux%M = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
