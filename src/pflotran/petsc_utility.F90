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

  interface PUMSetValues
    module procedure PUMSetValues1
    module procedure PUMSetValues2
    module procedure PUMSetValues3
  end interface

  interface PUMSetValuesLocal
    module procedure PUMSetValuesLocal1
    module procedure PUMSetValuesLocal2
    module procedure PUMSetValuesLocal3
  end interface

  interface PUCast
    module procedure PUCastBoolean
  end interface

  public :: PetUtilMatSVBL, &
            PetUtilVecSVBL, &
            PetUtilLoadVec, &
            PetUtilUnloadVec, &
            PetscUtilCompareMatrices

  public :: PUMSetValue, &
            PUMSetValues, &
            PUMSetValuesLocal, &
            PUMSetValuesBlocked, &
            PUMSetValuesBlockedLocal

  public :: PetscTestFile, &
            PUCast

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

  call MatSetValuesBlockedLocal(A,1,[irow-1],1,[icol-1], &
                                reshape(ndof_mat,(/size(ndof_mat)/)), &
                                ADD_VALUES,ierr);CHKERRQ(ierr)

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

  call VecGetArray(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = iarray(:)
  call VecRestoreArray(vec,vec_ptr,ierr);CHKERRQ(ierr)

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

  call VecGetArray(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = rarray(:)
  call VecRestoreArray(vec,vec_ptr,ierr);CHKERRQ(ierr)

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

  call VecGetArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)
  iarray(:) = nint(vec_ptr(:))
  call VecRestoreArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)

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

  call VecGetArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)
  rarray(:) = vec_ptr(:)
  call VecRestoreArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilUnloadVecReal

! ************************************************************************** !

subroutine PetscUtilCompareMatrices(A,B,nL2G,nG2A,rtol,row_rtol,option)
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
  PetscReal :: rtol
  PetscReal :: row_rtol
  type(option_type) :: option

  PetscInt :: nrow, ncol
  PetscInt, pointer :: irow_ptr(:), icol_ptr(:)
  PetscReal :: values_a(1000), values_b(1000)
  PetscReal :: value_a, value_b, scale, row_scale
  PetscBool :: success
  PetscInt :: i, j
  PetscInt :: irow(1)
  PetscInt :: irow_local0, icol_local0
  PetscInt :: irow_cell_local, icol_cell_local
  PetscInt :: row_cell_natural, col_cell_natural
  PetscInt :: irow_natural, icol_natural
  PetscInt :: block_size
  PetscErrorCode :: ierr

  call MatGetRowIJ(A,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                   nrow,irow_ptr,icol_ptr,success,ierr);CHKERRQ(ierr)
  if (.not.success) then
    option%io_buffer = 'Error returned from MatGetRowIJF90() in &
      &PetscUtilCompareMatrices. I believe that the type of matrix A is not &
      &supported by that routine.'
    call PrintErrMsg(option)
  endif
  call MatGetBlockSize(A,block_size,ierr);CHKERRQ(ierr)

  do i = 1, nrow
    irow(1) = i-1
    irow_local0 = i-1
    irow_cell_local = irow_local0/block_size+1
    row_cell_natural = nG2A(nL2G(irow_cell_local))
    irow_natural = (row_cell_natural-1)*block_size+mod(i-1,block_size)+1
    ncol = irow_ptr(i+1)-irow_ptr(i)
    call MatGetValuesLocal(A,ONE_INTEGER,[irow],ncol, &
                           icol_ptr(irow_ptr(i)+1:irow_ptr(i+1)), &
                           values_a,ierr);CHKERRQ(ierr)
    call MatGetValuesLocal(B,ONE_INTEGER,[irow],ncol, &
                           icol_ptr(irow_ptr(i)+1:irow_ptr(i+1)), &
                           values_b,ierr);CHKERRQ(ierr)
    do j = 1, ncol
      value_a = values_a(j)
      value_b = values_b(j)
      row_scale = max(maxval(dabs(values_a(1:ncol))), &
                      maxval(dabs(values_b(1:ncol))))
      scale = max(dabs(value_a),dabs(value_b))
      if (scale > 0.d0) then
        if (dabs(value_a - value_b)/scale > rtol .and. &
            dabs(value_a - value_b)/row_scale > row_rtol) then
          icol_local0 = (icol_ptr(irow_ptr(i)+j)+1)-1
          icol_cell_local = icol_local0/block_size+1
          col_cell_natural = nG2A(nL2G(icol_cell_local))
          icol_natural = (col_cell_natural-1)*block_size+mod(j-1,block_size)+1
          option%io_buffer = 'Large difference - (' // &
            StringWrite(irow_natural) // ',' // &
            StringWrite(icol_natural) // ') : ' // &
            StringWrite(value_a) // ' vs ' // StringWrite(value_b)
          call PrintMsgByRank(option)
        endif
      endif
    enddo
  enddo
  call MatRestoreRowIJ(A,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                       nrow,irow_ptr,icol_ptr,success,ierr);CHKERRQ(ierr)

end subroutine PetscUtilCompareMatrices

! ************************************************************************** !

subroutine PUMSetValue(A,irow,icol,scalar, &
                       insert_mode,ierr)
  !
  ! Maps to MatSetValue()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: irow
  PetscInt :: icol
  PetscReal :: scalar
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValue(A,irow,icol,scalar,insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValue

! ************************************************************************** !

subroutine PUMSetValues1(A,nrow,irow,ncol,icol,scalar, &
                         insert_mode,ierr)
  !
  ! Maps to MatSetValues()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: scalar
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValues(A,nrow,[irow],ncol,[icol],[scalar], &
                    insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValues1

! ************************************************************************** !

subroutine PUMSetValues2(A,nrow,irow,ncol,icol,array1d, &
                         insert_mode,ierr)
  !
  ! Maps to MatSetValues()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: array1d(:)
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValues(A,nrow,[irow],ncol,icol,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValues2

! ************************************************************************** !

subroutine PUMSetValues3(A,nrow,irow,ncol,icol,array1d, &
                         insert_mode,ierr)
  !
  ! Maps to MatSetValues()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow(:)
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: array1d(:)
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValues(A,nrow,irow,ncol,icol,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValues3

! ************************************************************************** !

subroutine PUMSetValuesLocal1(A,nrow,irow,ncol,icol,scalar, &
                              insert_mode,ierr)
  !
  ! Maps to MatSetValuesLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: scalar
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValuesLocal(A,nrow,[irow],ncol,[icol],[scalar], &
                         insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValuesLocal1

! ************************************************************************** !

subroutine PUMSetValuesLocal2(A,nrow,irow,ncol,icol,array2d, &
                              insert_mode,ierr)
  !
  ! Maps to MatSetValuesLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: array2d(:,:)
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValuesLocal(A,nrow,[irow],ncol,[icol], &
                         reshape(array2d,(/size(array2d)/)), &
                         insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValuesLocal2

! ************************************************************************** !

subroutine PUMSetValuesLocal3(A,nrow,irow,ncol,icol,array1d, &
                              insert_mode,ierr)
  !
  ! Maps to MatSetValuesLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: array1d(:)
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValuesLocal(A,nrow,[irow],ncol,icol,array1d, &
                         insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValuesLocal3

! ************************************************************************** !

subroutine PUMSetValuesBlocked(A,nrow,irow,ncol,icol,matrix_block, &
                                    insert_mode,ierr)
  !
  ! Maps to MatSetValuesBlocked()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: matrix_block(:,:)
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValuesBlocked(A,nrow,[irow],ncol,[icol], &
                           reshape(matrix_block,(/size(matrix_block)/)), &
                           insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValuesBlocked

! ************************************************************************** !

subroutine PUMSetValuesBlockedLocal(A,nrow,irow,ncol,icol,matrix_block, &
                                    insert_mode,ierr)
  !
  ! Maps to MatSetValuesBlockedLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: matrix_block(:,:)
  type(einsertmode) :: insert_mode
  PetscErrorCode :: ierr

  call MatSetValuesBlockedLocal(A,nrow,[irow],ncol,[icol], &
                                reshape(matrix_block,(/size(matrix_block)/)), &
                                insert_mode,ierr);CHKERRQ(ierr)

end subroutine PUMSetValuesBlockedLocal

! ************************************************************************** !

subroutine PetscTestFile(filename,char,flag,ierr)

  implicit none

  character(len=*) :: filename
  character(len=*) :: char
  PetscBool :: flag
  PetscErrorCode :: ierr

  inquire(file=trim(filename), exist=flag)
  ierr = 0

end subroutine PetscTestFile

! ************************************************************************** !

function PUCastBoolean(flag)

  implicit none

  logical :: flag

  PetscBool :: PUCastBoolean

  PUCastBoolean = flag

end function PUCastBoolean

end module Petsc_Utility_module
