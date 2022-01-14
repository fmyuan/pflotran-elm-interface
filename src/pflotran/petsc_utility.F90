module Petsc_Utility_module

 ! Wrapper routines for simplifying the interface between PFLOTRAN and 
 ! PETSc

#include "petsc/finclude/petscmat.h"
  use petscmat
 
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: PetUtilMatSVBL, &
            PetUtilVecSVBL
 
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

  call MatSetValuesBlockedLocal(A,1,irow-1,1,icol-1,matrix_block,ADD_VALUES, &
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
  PetscReal :: array(ndof)
  PetscInt :: icell
  PetscReal :: array_block(:)
  PetscBool :: replace

  PetscInt :: istart, iend
  PetscErrorCode :: ierr

  iend = icell*ndof
  istart = icell-ndof+1
  if (replace) then
    array(istart:iend) = array_block(1:ndof)
  else
    array(istart:iend) = array(istart:iend) + array_block(1:ndof)
  endif

end subroutine PetUtilVecSVBL

end module Petsc_Utility_module
