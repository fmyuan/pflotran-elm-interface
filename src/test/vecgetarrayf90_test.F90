program test
#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  Vec :: vec
  
  PetscReal, pointer :: vec_ptr(:)
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  call VecCreate(PETSC_COMM_SELF,vec,ierr);CHKERRQ(ierr)
  if (rank == 0) then
    call VecSetSizes(vec,16,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  else
    call VecSetSizes(vec,0,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  endif
  call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  if (rank == 0) then
    print *, vec_ptr(1)
  endif
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(vec, ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr);CHKERRQ(ierr)
 
end program test
