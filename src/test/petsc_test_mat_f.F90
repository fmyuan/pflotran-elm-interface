program test
#include "petsc/finclude/petscmat.h"
  use petscmat
!  use Create_Matrix
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  Mat :: A
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test mat program'

  call CreateMatrix(A,4)
  call MatDestroy(A,ierr);CHKERRA(ierr)
  if (rank == 0) print *
  if (rank == 0) print *, 'This line should not be reached due to the error &
    &above.'
  if (rank == 0) print *

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr)
 
end program test

subroutine CreateMatrix(A,n)
  use petscmat
  Mat :: A
  PetscInt :: n
  PetscReal :: value_
  call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, &
                    n,n, &
                    1,PETSC_NULL_INTEGER, &
                    2,PETSC_NULL_INTEGER,A,ierr);CHKERRQ(ierr)
  call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  value_ = 1.
  call MatSetValue(A,0,0,value_,INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A,0,1,value_,INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValue(A,0,2,value_,INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
end subroutine CreateMatrix
