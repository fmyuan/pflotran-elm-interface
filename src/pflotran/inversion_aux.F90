module Inversion_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module

  implicit none

  private

#if 0
  type, public :: inversion_auxvar_type
  contains
  end type inversion_auxvar_type
#endif

  type, public :: inversion_aux_type
    Mat :: JsensitivityT
    PetscInt :: num_aux
!    class(inversion_auxvar_type), pointer :: auxvars(:)
    PetscReal, pointer :: dFluxdIntConn(:,:)
    PetscReal, pointer :: dFluxdBCConn(:,:)
    PetscInt, pointer :: cell_to_internal_connection(:,:)
    PetscInt, pointer :: cell_to_bc_connection(:,:)
    Mat, pointer :: dMdK(:)
    Vec, pointer :: dbdK(:)
  end type inversion_aux_type

  public :: InversionAuxCreate, &
!            InversionAuxVarInit, &
!            InversionAuxVarStrip, &
            InversionAuxDestroy

contains

! ************************************************************************** !

function InversionAuxCreate()
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !

  use Option_module

  implicit none

  type(inversion_aux_type), pointer :: InversionAuxCreate

  type(inversion_aux_type), pointer :: aux

  allocate(aux)
!  nullify(aux%auxvars)

  aux%num_aux = 0
  nullify(aux%dFluxdIntConn)
  nullify(aux%cell_to_internal_connection)
  nullify(aux%cell_to_bc_connection)
  nullify(aux%dMdK)
  nullify(aux%dbdK)

  aux%JsensitivityT = PETSC_NULL_MAT

  InversionAuxCreate => aux

end function InversionAuxCreate

#if 0
! ************************************************************************** !

subroutine InversionAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Option_module

  implicit none

  class(inversion_auxvar_type) :: auxvar
  type(option_type) :: option

end subroutine InversionAuxVarInit

! ************************************************************************** !

subroutine InversionAuxVarStrip(auxvar)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_auxvar_type) :: auxvar

end subroutine InversionAuxVarStrip
#endif

! ************************************************************************** !

subroutine InversionAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(inversion_aux_type), pointer :: aux

  PetscInt :: iaux
  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

#if 0
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call InversionAuxVarStrip(aux%auxvars(iaux))
    enddo
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
#endif
  call DeallocateArray(aux%dFluxdIntConn)
  call DeallocateArray(aux%dFluxdBCConn)
  call DeallocateArray(aux%cell_to_internal_connection)
  call DeallocateArray(aux%cell_to_bc_connection)
  if (associated(aux%dMdK)) then
    do iaux = 1, size(aux%dMdK)
      call MatDestroy(aux%dMdK(iaux),ierr);CHKERRQ(ierr)
      call VecDestroy(aux%dbdK(iaux),ierr);CHKERRQ(ierr)
    enddo
    deallocate(aux%dMdK)
    nullify(aux%dMdK)
    deallocate(aux%dbdK)
    nullify(aux%dbdK)
  endif

  ! these objects are destroyed elsewhere, do not destroy
  aux%JsensitivityT = PETSC_NULL_MAT

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
