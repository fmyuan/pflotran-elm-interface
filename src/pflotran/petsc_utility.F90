module Petsc_Utility_module

 ! Wrapper routines for simplifying the interface between PFLOTRAN and
 ! PETSc

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module

  implicit none

  private

  interface PetUtilLoadVec
    module procedure PetUtilLoadVecInt
    module procedure PetUtilLoadVecReal
  end interface

  interface PetUtilUnloadVec
    module procedure PetUtilUnloadVecInt
    module procedure PetUtilUnloadVecReal
  end interface

  public :: PetUtilMatSVBL, &
            PetUtilVecSVBL, &
            PetUtilLoadVec, &
            PetUtilUnloadVec, &
            PetscUtilCompareMatrices

contains

! ************************************************************************** !

subroutine PetUtilMatSVBL(A,irow,icol,matrix_block,ndof)
  !
  ! Maps values from ndof*ndof sized array into ncell*ndof blocked matrix
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/22
  !
  implicit none

  Mat :: A
  PetscInt :: irow
  PetscInt :: icol
  PetscReal :: matrix_block(:,:)
  PetscInt :: ndof

  PetscReal :: ndof_mat(ndof,ndof)
  PetscErrorCode :: ierr

  ndof_mat = matrix_block(1:ndof,1:ndof)

  call MatSetValuesBlockedLocal(A,1,irow-1,1,icol-1,ndof_mat,ADD_VALUES, &
                                ierr);CHKERRQ(ierr)

end subroutine PetUtilMatSVBL

! ************************************************************************** !

subroutine PetUtilVecSVBL(array,icell,array_block,ndof,replace)
  !
  ! Maps values from ndof sized array into ncell*ndof blocked array
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/22
  !
  implicit none

  PetscInt :: ndof
  PetscReal :: array(:)
  PetscInt :: icell
  PetscReal :: array_block(:)
  PetscBool :: replace

  PetscInt :: istart, iend

  iend = icell*ndof
  istart = iend-ndof+1
  if (replace) then
    array(istart:iend) = array_block(1:ndof)
  else
    array(istart:iend) = array(istart:iend) + array_block(1:ndof)
  endif

end subroutine PetUtilVecSVBL

! ************************************************************************** !

subroutine PetUtilLoadVecInt(vec,iarray)
  !
  ! Loads values from array into a vec (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscInt :: iarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = iarray(:)
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilLoadVecInt

! ************************************************************************** !

subroutine PetUtilLoadVecReal(vec,rarray)
  !
  ! Loads values from array into a vec (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscReal :: rarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = rarray(:)
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilLoadVecReal

! ************************************************************************** !

subroutine PetUtilUnloadVecInt(vec,iarray)
  !
  ! Unloads values from a vec into array (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscInt :: iarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayReadF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  iarray(:) = int(vec_ptr(:)+1.d-5)
  call VecRestoreArrayReadF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilUnloadVecInt

! ************************************************************************** !

subroutine PetUtilUnloadVecReal(vec,rarray)
  !
  ! Unloads values from a vec into array (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscReal :: rarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayReadF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  rarray(:) = vec_ptr(:)
  call VecRestoreArrayReadF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilUnloadVecReal

! ************************************************************************** !

subroutine PetscUtilCompareMatrices(A,B,nL2G,nG2A,option)
  !
  ! Compares values in matrices pinpointing differences.
  !
  ! Author: Glenn Hammond
  ! Date: 01/02/25
  !
  use Option_module
  use String_module

  implicit none

  Mat :: A, B
  PetscInt :: nL2G(:) ! from grid%nL2G
  PetscInt :: nG2A(:)
  type(option_type) :: option

  PetscInt :: nrow, ncol
  PetscInt, pointer :: irow_ptr(:), icol_ptr(:)
  PetscReal :: values_a(1,1000), values_b(1,1000)
  PetscReal :: value_a, value_b, scale, row_scale
  PetscBool :: success
  PetscInt :: i, j
  PetscInt :: irow(1)
  PetscErrorCode :: ierr

  call MatGetRowIJF90(A,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                      nrow,irow_ptr,icol_ptr,success,ierr);CHKERRQ(ierr)
  do i = 1, nrow
    irow(1) = i-1
    ncol = irow_ptr(i+1)-irow_ptr(i)
    call MatGetValuesLocal(A,ONE_INTEGER,irow,ncol, &
                           icol_ptr(irow_ptr(i)+1:irow_ptr(i+1)), &
                           values_a,ierr);CHKERRQ(ierr)
    call MatGetValuesLocal(B,ONE_INTEGER,irow,ncol, &
                           icol_ptr(irow_ptr(i)+1:irow_ptr(i+1)), &
                           values_b,ierr);CHKERRQ(ierr)
    do j = 1, ncol
      value_a = values_a(1,j)
      value_b = values_b(1,j)
      row_scale = max(maxval(dabs(values_a(1,1:ncol))), &
                      maxval(dabs(values_b(1,1:ncol))))
      scale = max(dabs(value_a),dabs(value_b))
      if (scale > 0.d0) then
        if (dabs(value_a - value_b)/scale > 1.d-1 .and. &
            dabs(value_a - value_b)/row_scale > 1.d-8) then
          option%io_buffer = 'Large difference - (' // &
            StringWrite(nG2A(nL2G(i))) // ',' // &
            StringWrite(nG2A(icol_ptr(j)+1)) // ') : ' // &
            StringWrite(value_a) // ' vs ' // StringWrite(value_b)
          call PrintMsgByRank(option)
        endif
      endif
    enddo
  enddo
  call MatRestoreRowIJF90(A,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                          nrow,irow_ptr,icol_ptr,success,ierr);CHKERRQ(ierr)

end subroutine PetscUtilCompareMatrices

end module Petsc_Utility_module
