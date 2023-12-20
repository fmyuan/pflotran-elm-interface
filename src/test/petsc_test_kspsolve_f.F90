program test
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscdmda.h"
#include "petsc/finclude/petscksp.h"
  use petscsys
  use petscmat
  use petscvec
  use petscdmda
  use petscksp
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  PetscViewer :: viewer
  
  character(len=32) :: filename
  
  PetscInt :: ndof
  PetscInt :: ncell

  Mat :: A
  Vec :: x
  Vec :: b

  DM :: da

  KSP :: ksp

  PetscInt :: ncell_local
  PetscInt :: offset
  PetscInt :: dummy
  PetscInt :: i
  PetscReal :: value
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  ncell = 3
  ndof = 11

!  call DMDACreate1D(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,ncell,ndof,1, &
!                    PETSC_NULL_INTEGER,da,ierr);CHKERRQ(ierr)
  call DMDACreate3D(PETSC_COMM_WORLD,DM_BOUNDARY_NONE, &
                    DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR, &
                    ncell,1,1, &
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
                    ndof,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    da,ierr);CHKERRQ(ierr)
  call DMSetup(da,ierr);CHKERRQ(ierr)
  call DMDAGetCorners(da,offset,dummy,dummy,ncell_local,dummy,dummy,ierr);CHKERRQ(ierr)
  call DMCreateGlobalVector(da,x,ierr);CHKERRQ(ierr)
  call VecDuplicate(x,b,ierr);CHKERRQ(ierr)
  call DMSetMatType(da,MATBAIJ,ierr);CHKERRQ(ierr)
  call DMCreateMatrix(da,A,ierr);CHKERRQ(ierr)

  call VecZeroEntries(b,ierr);CHKERRQ(ierr)
  if (rank == 0) then
    value = 6.66134d-16
    call VecSetValue(b,0,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
    value = -value
    call VecSetValue(b,ndof,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
  endif
  call VecAssemblyBegin(b,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(b,ierr);CHKERRQ(ierr)

  do i = offset*ndof, (offset+ncell_local)*ndof-1
    value = 1333.33
    call MatSetValue(A,i,i,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
    if (i > ndof-1) then
      value = -7.5d-7
      call MatSetValue(A,i,i-ndof,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
    endif
    if (i < (ncell-1)*ndof) then
      value = -7.5d-7
      call MatSetValue(A,i,i+ndof,value,INSERT_VALUES,ierr);CHKERRQ(ierr)
    endif
  enddo
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  filename = 'Af.txt'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  filename = 'bf.txt'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr);CHKERRQ(ierr)
  call VecView(b,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr);CHKERRQ(ierr)
  call KSPSetType(ksp,KSPBCGS,ierr);CHKERRQ(ierr)
  call KSPSetFromOptions(ksp,ierr);CHKERRQ(ierr)
  call KSPSetOperators(ksp,A,A,ierr);CHKERRQ(ierr)
  call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)

  filename = 'xf.txt'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr);CHKERRQ(ierr)
  call VecView(x,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call VecDestroy(x,ierr);CHKERRQ(ierr)
  call VecDestroy(b,ierr);CHKERRQ(ierr)
  call MatDestroy(A,ierr);CHKERRQ(ierr)
  call DMDestroy(da,ierr);CHKERRQ(ierr)
  call KSPDestroy(ksp,ierr);CHKERRQ(ierr)
  call PetscFinalize(ierr);CHKERRQ(ierr)
 
  if (rank == 0) print *, 'End of Fortran90 test program'
 
end program test
