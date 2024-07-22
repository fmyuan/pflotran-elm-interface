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
            PetUtilUnloadVec

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

end module Petsc_Utility_module
