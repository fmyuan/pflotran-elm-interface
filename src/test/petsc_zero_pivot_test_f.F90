program test
#include "petsc/finclude/petscksp.h"
  use petscksp
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  Vec :: x
  Vec :: b
  Mat :: A
  KSP :: ksp
  PC :: pc
  PetscReal :: value
  PetscReal :: tolerance
  PetscInt :: i, n
  PetscInt :: offset
  PetscReal, pointer :: x_ptr(:)
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  n = 2/size
  call MatCreateAIJ(PETSC_COMM_WORLD,n,n, &
                    PETSC_DETERMINE,PETSC_DETERMINE, &
                    1,PETSC_NULL_INTEGER, &
                    0,PETSC_NULL_INTEGER,A,ierr);CHKERRQ(ierr)
  call VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DETERMINE,x,ierr);CHKERRQ(ierr)
  call VecDuplicate(x,b,ierr);CHKERRQ(ierr)
  offset = rank * n
  value = 1.d-16
  do i = 1, n
    call MatSetValue(A,i-1+offset,i-1+offset,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
    call VecSetValue(b,i-1+offset,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr);CHKERRQ(ierr)
  call KSPSetErrorIfNotConverged(ksp,PETSC_TRUE,ierr);CHKERRQ(ierr)
  call KSPSetOperators(ksp,A,A,ierr);CHKERRQ(ierr)
  call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
  call KSPSetFromOptions(ksp,ierr);CHKERRQ(ierr)
  call PCSetFromOptions(pc,ierr);CHKERRQ(ierr)
!  call KSPSetup(ksp,ierr);CHKERRQ(ierr)

  call PCFactorSetShiftType(pc,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
  tolerance = 1.d-20
  call PCFactorSetZeroPivot(pc,tolerance,ierr);CHKERRQ(ierr)

!  call KSPSetup(ksp,ierr);CHKERRQ(ierr)
  call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(x,x_ptr,ierr);CHKERRQ(ierr)
  print *, 'These values should be ~1: ', x_ptr 
  call VecRestoreArrayF90(x,x_ptr,ierr);CHKERRQ(ierr)

  call KSPDestroy(ksp,ierr);CHKERRQ(ierr)
  call MatDestroy(A,ierr);CHKERRQ(ierr)
  call VecDestroy(x,ierr);CHKERRQ(ierr)
  call VecDestroy(b,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr);CHKERRQ(ierr)
 
end program test
