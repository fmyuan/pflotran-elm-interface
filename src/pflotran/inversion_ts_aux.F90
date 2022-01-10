module Inversion_TS_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_ts_aux_type
    PetscInt :: timestep
    PetscReal :: time
    Mat :: dResdu             ! copy of Jacobian: df(u)/du
    Mat :: dResdparam         ! matrix sotring df(u)/dparameter
    Vec, pointer :: lambda(:) ! arrays storing dg(u)/df(u)
    ! derivative of residual wrt unknown at old time level (k)
    PetscReal, pointer :: dRes_du_k(:)  ! array storing df(u)^k+1/du^k
    PetscReal, pointer :: dFluxdIntConn(:,:)
    PetscReal, pointer :: dFluxdBCConn(:,:)
    Vec :: solution
    type(inversion_mat_vec_pointer_type), pointer :: mat_vec_solution_ptr
    type(inversion_ts_aux_type), pointer :: prev
    type(inversion_ts_aux_type), pointer :: next
  end type inversion_ts_aux_type

  type, public :: inversion_mat_vec_pointer_type
    Mat :: M
    Vec :: solution
  end type inversion_mat_vec_pointer_type

  public :: InversionTSAuxCreate, &
            InvTSAuxAllocate, &
            InvTSAuxAllocateFluxCoefArrays, &
            InvTSAuxStoreCopyGlobalMatVecs, &
            InversionTSAuxListDestroy, &
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
  aux%dResdu = PETSC_NULL_MAT
  aux%dResdparam = PETSC_NULL_MAT
  aux%solution = PETSC_NULL_VEC
  nullify(aux%lambda)

  nullify(aux%mat_vec_solution_ptr)
  nullify(aux%dRes_du_k)
  nullify(aux%dFluxdIntConn)
  nullify(aux%dFluxdBCConn)
  nullify(aux%prev)
  nullify(aux%next)

  if (associated(prev_ts_aux)) then
    aux%timestep = prev_ts_aux%timestep + 1
    aux%mat_vec_solution_ptr => prev_ts_aux%mat_vec_solution_ptr
    prev_ts_aux%next => aux
    aux%prev => prev_ts_aux
    call InvTSAuxAllocate(aux,size(prev_ts_aux%dRes_du_k))
    if (associated(prev_ts_aux%dFluxdIntConn)) then
      call InvTSAuxAllocateFluxCoefArrays(aux, &
                                          size(prev_ts_aux%dFluxdIntConn,2), &
                                          size(prev_ts_aux%dFluxdBCConn,2))
    endif
  else
    aux%timestep = 1
    allocate(aux%mat_vec_solution_ptr)
    aux%mat_vec_solution_ptr%M = PETSC_NULL_MAT
    aux%mat_vec_solution_ptr%solution = PETSC_NULL_VEC
  endif

  InversionTSAuxCreate => aux

end function InversionTSAuxCreate

! ************************************************************************** !

subroutine InvTSAuxAllocate(aux,num_unknown)
  !
  ! Allocated array holding Jacobian and dResdk matrix coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_ts_aux_type), pointer :: aux
  PetscInt :: num_unknown
  PetscErrorCode :: ierr

  call MatDuplicate(aux%mat_vec_solution_ptr%M,MAT_SHARE_NONZERO_PATTERN, &
                    aux%dResdparam,ierr);CHKERRQ(ierr)
  allocate(aux%dRes_du_k(num_unknown))
  aux%dRes_du_k = 0.d0

end subroutine InvTSAuxAllocate

! ************************************************************************** !

subroutine InvTSAuxAllocateFluxCoefArrays(aux,num_internal,num_boundary)
  !
  ! Allocated array holding flux coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_ts_aux_type), pointer :: aux
  PetscInt :: num_internal
  PetscInt :: num_boundary
  PetscErrorCode :: ierr

  allocate(aux%dFluxdIntConn(6,num_internal))
  aux%dFluxdIntConn = 0.d0
  allocate(aux%dFluxdBCConn(2,num_boundary))
  aux%dFluxdBCConn = 0.d0

end subroutine InvTSAuxAllocateFluxCoefArrays

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
                    aux%dResdu,ierr);CHKERRQ(ierr)
  if (aux%mat_vec_solution_ptr%solution /= PETSC_NULL_VEC) then
    call VecDuplicate(aux%mat_vec_solution_ptr%solution, &
                      aux%solution,ierr);CHKERRQ(ierr)
    call VecCopy(aux%mat_vec_solution_ptr%solution,aux%solution, &
                ierr);CHKERRQ(ierr)
  endif

end subroutine InvTSAuxStoreCopyGlobalMatVecs

! ************************************************************************** !

subroutine InversionTSAuxListDestroy(inversion_ts_aux_list,print_msg)
  !
  ! Deallocates a inversion auxiliary timestep list object
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/21
  !
  type(inversion_ts_aux_type), pointer :: inversion_ts_aux_list
  PetscBool :: print_msg

  type(inversion_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_ts_aux_type), pointer :: next_inversion_ts_aux

  cur_inversion_ts_aux => inversion_ts_aux_list
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    next_inversion_ts_aux => cur_inversion_ts_aux%next
#if 0
    if (.not.associated(next_inversion_ts_aux)) then
      if (print_msg) then
        print *
        print *, '  Last Inversion TS Aux timestep and time: ', &
          cur_inversion_ts_aux%timestep, cur_inversion_ts_aux%time
        print *
      endif
    endif
#endif
    call InversionTSAuxDestroy(cur_inversion_ts_aux)
    cur_inversion_ts_aux => next_inversion_ts_aux
  enddo

  nullify(inversion_ts_aux_list)

end subroutine InversionTSAuxListDestroy

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
  if (aux%dResdu /= PETSC_NULL_MAT) then
    call MatDestroy(aux%dResdu,ierr);CHKERRQ(ierr)
  endif
  if (aux%dResdparam /= PETSC_NULL_MAT) then
    call MatDestroy(aux%dResdparam,ierr);CHKERRQ(ierr)
  endif
  if (aux%solution /= PETSC_NULL_VEC) then
    call VecDestroy(aux%solution,ierr);CHKERRQ(ierr)
  endif
  if (associated(aux%lambda)) then
    call VecDestroyVecs(size(aux%lambda),aux%lambda,ierr);CHKERRQ(ierr)
    nullify(aux%lambda)
  endif
  call DeallocateArray(aux%dRes_du_k)
  call DeallocateArray(aux%dFluxdIntConn)
  call DeallocateArray(aux%dFluxdBCConn)

  deallocate(aux)
  nullify(aux)

end subroutine InversionTSAuxDestroy

end module Inversion_TS_Aux_module
