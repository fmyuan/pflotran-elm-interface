program test

  use petscmat

  implicit none

#include "petsc/finclude/petscmat.h"

  PetscMPIInt :: size
  PetscMPIInt :: rank
  Mat :: A
  PetscInt, parameter :: nx = 4
  PetscInt, parameter :: ny = 4
  PetscInt, parameter :: n_local = 8
  PetscInt :: i, j, ii
  PetscInt :: d_nnz(8)
  PetscInt :: o_nnz(8)
  PetscInt :: j_offset
  PetscErrorCode :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test mat program'

  call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATBAIJ,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,n_local,n_local,PETSC_DETERMINE,PETSC_DETERMINE, &
                   ierr);CHKERRQ(ierr)
  o_nnz = 0
  d_nnz = 0
  if (rank == 0) then
    j_offset = 0
    d_nnz(1) = 3; o_nnz(1) = 0
    d_nnz(2) = 4; o_nnz(2) = 0
    d_nnz(3) = 4; o_nnz(3) = 0
    d_nnz(4) = 3; o_nnz(4) = 0
    d_nnz(5) = 3; o_nnz(5) = 1
    d_nnz(6) = 4; o_nnz(6) = 1
    d_nnz(7) = 4; o_nnz(7) = 1
    d_nnz(8) = 3; o_nnz(8) = 1
  else
    j_offset = 2
    d_nnz(1) = 3; o_nnz(5) = 1
    d_nnz(2) = 4; o_nnz(6) = 1
    d_nnz(3) = 4; o_nnz(7) = 1
    d_nnz(4) = 3; o_nnz(8) = 1
    d_nnz(5) = 3; o_nnz(1) = 0
    d_nnz(6) = 4; o_nnz(2) = 0
    d_nnz(7) = 4; o_nnz(3) = 0
    d_nnz(8) = 3; o_nnz(4) = 0
  endif
  call MatMPIBAIJSetPreallocation(A,1,PETSC_NULL_INTEGER,d_nnz, &
                                  PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
!  call MatSetUp(A,ierr);CHKERRQ(ierr)
  
if (rank == 1) then
  ii = 9
  call SetValue(A,ii,ii-nx)
endif

if (rank == 3) then
  do j = 1+j_offset,2+j_offset
    do i = 1, nx
      ii = i + (j-1)*nx
      if (j > 1) then
        call SetValue(A,ii,ii-nx)
      endif
      if (i > 1) then
        call SetValue(A,ii,ii-1)
      endif
      call SetValue(A,ii,ii)
      if (i < nx) then
        call SetValue(A,ii,ii+1)
      endif
      if (j < ny) then
        call SetValue(A,ii,ii+nx)
      endif
    enddo
  enddo
endif

print *, rank, 'here0'
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
print *, rank, 'here1'

#if 0
  if (rank == 0) then
    call SetValue(A,14,2)
  endif
#endif

  call MatDestroy(A,ierr);CHKERRA(ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr)
 
end program test

subroutine SetValue(A,irow,icol)

  use petscmat

  implicit none

#include "petsc/finclude/petscmat.h"

  Mat :: A
  PetscInt :: irow
  PetscInt :: icol

  PetscInt :: i
  PetscInt :: j
  PetscReal :: values(1,1)
  PetscErrorCode :: ierr

  values = 1.d0
  i = irow-1
  j = icol-1

  print *, irow, icol
  call MatSetValuesBlocked(A,1,i,1,j,values,ADD_VALUES,ierr);CHKERRQ(ierr)

end subroutine SetValue
