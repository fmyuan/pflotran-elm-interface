program test
#include "finclude/petscksp.h"
  use petscksp
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  Vec :: b
  Vec :: x
  Mat :: A
  KSP :: ksp
  PetscViewer :: viewer
  
  character(len=32) :: filename
  
  PetscInt :: ndof
  PetscInt :: my_ndof
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  ndof = 3
  call VecCreateSeq(PETSC_COMM_WORLD, ndof, b, ierr);CHKERRQ(ierr)
  call VecDuplicate(b, x, ierr);CHKERRQ(ierr)
  call MatCreateSeqAIJ(PETSC_COMM_WORLD, ndof, ndof, ndof, PETSC_NULL_INTEGER, &
                       A, ierr);CHKERRQ(ierr)

  call VecZeroEntries(b, ierr);CHKERRQ(ierr)
  call VecZeroEntries(x, ierr);CHKERRQ(ierr)
  call VecSetValue(b, 0, -0.00555083d0, INSERT_VALUES, ierr);CHKERRQ(ierr)
  call VecSetValue(b, 1, 0.d0, INSERT_VALUES, ierr);CHKERRQ(ierr)
  call VecSetValue(b, 2, 0.00555083d0, INSERT_VALUES, ierr);CHKERRQ(ierr)
  call MatZeroEntries(A, ierr);CHKERRQ(ierr)
  call MatSetValue(A, 0, 0, 5.55083d-7, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A, 0, 1, -5.55083d-7, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A, 1, 0, -5.55083d-7, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A, 1, 1, 1.110166d-6, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A, 1, 2, 5.55083d-7, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A, 2, 1, -5.55083d-7, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A, 2, 2, 5.55083d-7, INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  print *, 'b:'
  call VecView(b, PETSC_VIEWER_STDOUT_SELF, ierr);CHKERRQ(ierr)
  print *, 'A:'
  call MatView(A, PETSC_VIEWER_STDOUT_SELF, ierr);CHKERRQ(ierr)

  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr);CHKERRQ(ierr)
  call KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN, ierr);CHKERRQ(ierr)
  call KSPSolve(ksp, b, x, ierr);CHKERRQ(ierr)

  print *, 'x:'
  call VecView(x, PETSC_VIEWER_STDOUT_SELF, ierr);CHKERRQ(ierr)

  call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)
  call KSPDestroy(ksp, ierr);CHKERRQ(ierr)
  call VecDestroy(b, ierr);CHKERRQ(ierr)
  call VecDestroy(x, ierr);CHKERRQ(ierr)
  call MatDestroy(A, ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr);CHKERRQ(ierr)
 
end program test
