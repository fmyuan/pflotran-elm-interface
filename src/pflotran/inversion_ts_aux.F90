module Inversion_TS_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Inversion_Coupled_Aux_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module

  implicit none

  private

  type, public :: inversion_forward_ts_aux_type
    PetscInt :: timestep
    PetscReal :: time
    Mat :: dResdu             ! copy of Jacobian: df(u)/du
    Mat :: dResdparam         ! matrix storing df(u)/dparameter
    ! derivative of residual wrt unknown at old time level (k)
    PetscReal, pointer :: dRes_du_k(:,:,:)  ! array storing df(u)^k+1/du^k
    type(inversion_forward_ts_aux_type), pointer :: prev
    type(inversion_forward_ts_aux_type), pointer :: next
  end type inversion_forward_ts_aux_type

  interface InversionForwardTSAuxCreate
    module procedure :: InversionForwardTSAuxCreate1
    module procedure :: InversionForwardTSAuxCreate2
  end interface InversionForwardTSAuxCreate


  public :: InversionForwardTSAuxCreate, &
            InvForwardTSAuxDestroyList, &
            InvForTSAuxDupForwardJacobian, &
            InvForwardTSAuxDestroyMatrices

contains

! ************************************************************************** !

function InversionForwardTSAuxCreate1(M_ptr,ndof,ncell)
  !
  ! Allocate and initialize auxiliary inversion time step object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  implicit none

  Mat :: M_ptr
  PetscInt :: ndof
  PetscInt :: ncell

  type(inversion_forward_ts_aux_type), pointer :: InversionForwardTSAuxCreate1

  type(inversion_forward_ts_aux_type), pointer :: aux

  allocate(aux)
  call InvForwardTSAuxInit(aux)
  aux%timestep = 1
  call InvForwardTSAuxAllocateMatrices(aux,M_ptr,ndof,ncell)

  InversionForwardTSAuxCreate1 => aux

end function InversionForwardTSAuxCreate1

! ************************************************************************** !

function InversionForwardTSAuxCreate2(prev_ts_aux)
  !
  ! Allocate and initialize auxiliary inversion time step object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  implicit none

  type(inversion_forward_ts_aux_type), pointer :: prev_ts_aux

  type(inversion_forward_ts_aux_type), pointer :: InversionForwardTSAuxCreate2

  type(inversion_forward_ts_aux_type), pointer :: aux

  allocate(aux)
  call InvForwardTSAuxInit(aux)

  aux%timestep = prev_ts_aux%timestep + 1
  prev_ts_aux%next => aux
  aux%prev => prev_ts_aux
  call InvForwardTSAuxAllocateMatrices(aux,aux%prev%dResdparam, &
                                       size(prev_ts_aux%dRes_du_k,1), &
                                       size(prev_ts_aux%dRes_du_k,3))

  InversionForwardTSAuxCreate2 => aux

end function InversionForwardTSAuxCreate2

! ************************************************************************** !

subroutine InvForwardTSAuxInit(aux)
  !
  ! Initializes the contents of inversion_forward_ts_aux_type
  !
  ! Author: Glenn Hammond
  ! Date: 11/28/22

  implicit none

  type(inversion_forward_ts_aux_type) :: aux

  aux%timestep = UNINITIALIZED_INTEGER
  aux%time = UNINITIALIZED_DOUBLE
  aux%dResdu = PETSC_NULL_MAT
  aux%dResdparam = PETSC_NULL_MAT

  nullify(aux%dRes_du_k)
  nullify(aux%prev)
  nullify(aux%next)

end subroutine InvForwardTSAuxInit

! ************************************************************************** !

subroutine InvForwardTSAuxAllocateMatrices(aux,M_ptr,ndof,ncell)
  !
  ! Allocated array holding Jacobian and dResdk matrix coefficients
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21

  type(inversion_forward_ts_aux_type), pointer :: aux
  Mat :: M_ptr
  PetscInt :: ndof
  PetscInt :: ncell
  PetscErrorCode :: ierr

  call MatDuplicate(M_ptr,MAT_SHARE_NONZERO_PATTERN,aux%dResdparam, &
                    ierr);CHKERRQ(ierr)
  allocate(aux%dRes_du_k(ndof,ndof,ncell))
  aux%dRes_du_k = 0.d0

end subroutine InvForwardTSAuxAllocateMatrices

! ************************************************************************** !

subroutine InvForTSAuxDupForwardJacobian(M_ptr,ts_aux)
  !
  ! Copies Jacobian matrix and solution vector
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22

  Mat :: M_ptr
  type(inversion_forward_ts_aux_type), pointer :: ts_aux

  PetscErrorCode :: ierr

  call MatDuplicate(M_ptr,MAT_COPY_VALUES,ts_aux%dResdu, &
                    ierr);CHKERRQ(ierr)

end subroutine InvForTSAuxDupForwardJacobian

! ************************************************************************** !

subroutine InvForwardTSAuxDestroyList(list,print_msg)
  !
  ! Deallocates a inversion auxiliary timestep list object
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/21
  !
  type(inversion_forward_ts_aux_type), pointer :: list
  PetscBool :: print_msg

  type(inversion_forward_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_forward_ts_aux_type), pointer :: next_inversion_ts_aux

  cur_inversion_ts_aux => list
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
    call InvForwardTSAuxDestroyMatrices(cur_inversion_ts_aux)
    cur_inversion_ts_aux => next_inversion_ts_aux
  enddo

  nullify(list)

end subroutine InvForwardTSAuxDestroyList

! ************************************************************************** !

subroutine InvForwardTSAuxDestroyMatrices(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 11/19/21
  !
  use Utility_module, only : DeallocateArray

  type(inversion_forward_ts_aux_type), pointer :: aux

  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  nullify(aux%prev)
  nullify(aux%next)

  if (aux%dResdu /= PETSC_NULL_MAT) then
    call MatDestroy(aux%dResdu,ierr);CHKERRQ(ierr)
  endif
  if (aux%dResdparam /= PETSC_NULL_MAT) then
    call MatDestroy(aux%dResdparam,ierr);CHKERRQ(ierr)
  endif
  call DeallocateArray(aux%dRes_du_k)

  deallocate(aux)
  nullify(aux)

end subroutine InvForwardTSAuxDestroyMatrices

end module Inversion_TS_Aux_module
