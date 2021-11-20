module Inversion_TS_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_ts_aux_type
    PetscInt :: timestep
    PetscReal :: time
    Mat :: M
    Vec :: solution
    PetscReal, pointer :: dFluxdIntConn(:,:)
    PetscReal, pointer :: dFluxdBCConn(:,:)
    Mat, pointer :: dMdK(:)
    Vec, pointer :: dbdK(:)
    type(inversion_mat_vec_pointer_type), pointer :: mat_vec_solution_ptr
    type(inversion_ts_aux_type), pointer :: prev
    type(inversion_ts_aux_type), pointer :: next
  end type inversion_ts_aux_type

  type, public :: inversion_mat_vec_pointer_type
    Mat :: M
    Vec :: solution
  end type inversion_mat_vec_pointer_type

  public :: InversionTSAuxCreate, &
            InvTSAuxAllocateFluxCoefArrays, &
            InvTSAuxAllocateMatsAndVecs, &
            InvTSAuxStoreCopyGlobalMatVecs, &
            InversionTSAuxDestroy

contains

! ************************************************************************** !

function InversionTSAuxCreate(prev_ts_aux)
  !
  ! Allocate and initialize auxiliary inversion time step object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  implicit none

  type(inversion_ts_aux_type), pointer :: prev_ts_aux

  type(inversion_ts_aux_type), pointer :: InversionTSAuxCreate

  type(inversion_ts_aux_type), pointer :: aux

  allocate(aux)

  aux%time = UNINITIALIZED_DOUBLE
  aux%M = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC

  nullify(aux%mat_vec_solution_ptr)
  nullify(aux%dFluxdIntConn)
  nullify(aux%dFluxdBCConn)
  nullify(aux%dMdK)
  nullify(aux%dbdK)
  nullify(aux%prev)
  nullify(aux%next)

  if (associated(prev_ts_aux)) then
    aux%timestep = prev_ts_aux%timestep + 1
    aux%mat_vec_solution_ptr => prev_ts_aux%mat_vec_solution_ptr
    prev_ts_aux%next => aux
    aux%prev => prev_ts_aux
    call InvTSAuxAllocateFluxCoefArrays(aux, &
                                        size(prev_ts_aux%dFluxdIntConn,2), &
                                        size(prev_ts_aux%dFluxdBCConn,2))
    if (associated(prev_ts_aux%dMdK)) then
      call InvTSAuxAllocateMatsAndVecs(aux,size(aux%dMdK), &
                                       prev_ts_aux%dMdK(1), &
                                       prev_ts_aux%dbdK(1))
    endif
  else
    aux%timestep = 1
    allocate(aux%mat_vec_solution_ptr)
  endif

  InversionTSAuxCreate => aux

end function InversionTSAuxCreate

! ************************************************************************** !

subroutine InvTSAuxStoreCopyGlobalMatVecs(aux)
  !
  ! Allocated array holding internal flux coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_ts_aux_type), pointer :: aux

  PetscErrorCode :: ierr

  call MatDuplicate(aux%mat_vec_solution_ptr%M,MAT_COPY_VALUES, &
                    aux%M,ierr);CHKERRQ(ierr)
  call VecDuplicate(aux%mat_vec_solution_ptr%solution, &
                    aux%solution,ierr);CHKERRQ(ierr)
  call VecCopy(aux%mat_vec_solution_ptr%solution,aux%solution, &
               ierr);CHKERRQ(ierr)

end subroutine InvTSAuxStoreCopyGlobalMatVecs

! ************************************************************************** !

subroutine InvTSAuxAllocateFluxCoefArrays(aux,num_internal,num_boundary)
  !
  ! Allocated array holding internal flux coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_ts_aux_type), pointer :: aux
  PetscInt :: num_internal
  PetscInt :: num_boundary

  allocate(aux%dFluxdIntConn(6,num_internal))
  aux%dFluxdIntConn = 0.d0
  allocate(aux%dFluxdBCConn(2,num_boundary))
  aux%dFluxdBCConn = 0.d0

end subroutine InvTSAuxAllocateFluxCoefArrays

! ************************************************************************** !

subroutine InvTSAuxAllocateMatsAndVecs(aux,num_mats_and_vecs,mat,vec)
  !
  ! Allocated array holding internal flux coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_ts_aux_type), pointer :: aux
  PetscInt :: num_mats_and_vecs
  Mat :: mat
  Vec :: vec

  PetscInt :: i
  PetscErrorCode :: ierr

  allocate(aux%dMdK(num_mats_and_vecs))
  aux%dMdK = PETSC_NULL_MAT
  allocate(aux%dbdK(num_mats_and_vecs))
  aux%dbdK = PETSC_NULL_VEC
  do i = 1, num_mats_and_vecs
    call MatDuplicate(mat,MAT_SHARE_NONZERO_PATTERN, &
                      aux%dMdK(i),ierr);CHKERRQ(ierr)
  enddo
  do i = 1, num_mats_and_vecs
    call VecDuplicate(vec,aux%dbdK(i),ierr);CHKERRQ(ierr)
  enddo

end subroutine InvTSAuxAllocateMatsAndVecs

! ************************************************************************** !

subroutine InversionTSAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  use Utility_module, only : DeallocateArray

  type(inversion_ts_aux_type), pointer :: aux

  PetscInt :: iaux
  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  nullify(aux%prev)
  nullify(aux%next)

  ! only destroy once
  if (aux%timestep == 1) then
    deallocate(aux%mat_vec_solution_ptr)
  endif
  nullify(aux%mat_vec_solution_ptr)
  call MatDestroy(aux%M,ierr);CHKERRQ(ierr)
  call VecDestroy(aux%solution,ierr);CHKERRQ(ierr)
  call DeallocateArray(aux%dFluxdIntConn)
  call DeallocateArray(aux%dFluxdBCConn)
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

  deallocate(aux)
  nullify(aux)

end subroutine InversionTSAuxDestroy

end module Inversion_TS_Aux_module
