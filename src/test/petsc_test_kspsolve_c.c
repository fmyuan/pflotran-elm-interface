#include <stdio.h>

#include <petscsys.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscviewer.h>

int main(int argc,char* argv[]){

  PetscMPIInt rank;
  PetscMPIInt size;

  PetscInt ndof;
  PetscInt ncell;

  Mat A;
  Vec x;
  Vec b;

  PetscViewer viewer;

  DM da;
  KSP ksp;

  PetscInt ncell_local;
  PetscInt offset;
  PetscInt dummy;
  PetscInt i;
  PetscReal value;

  DMBoundaryType bt = DM_BOUNDARY_NONE;
  DMDAStencilType stype = DMDA_STENCIL_STAR;
  

  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,(char *)0);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) printf("Beginning of C test program\n");

  ncell = 3;
  ndof = 11;

  DMDACreate3d(PETSC_COMM_WORLD,bt,bt,bt,stype,
               ncell,1,1, 
               PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 
               ndof,1,NULL,NULL,NULL,&da);

  DMSetUp(da);
  DMDAGetCorners(da,&offset,&dummy,&dummy,&ncell_local,&dummy,&dummy);
  DMCreateGlobalVector(da,&x);
  VecDuplicate(x,&b);
  DMSetMatType(da,MATBAIJ);
  DMCreateMatrix(da,&A);

  for (i=offset*ndof;i<(offset+ncell_local)*ndof;i++) {
    value = 1333.33;
    MatSetValue(A,i,i,value,INSERT_VALUES);
    if (i > ndof-1) {
      value = -7.5e-7;
      MatSetValue(A,i,i-ndof,value,INSERT_VALUES);
    }
    if (i < (ncell-1)*ndof) {
      value = -7.5e-7;
      MatSetValue(A,i,i+ndof,value,INSERT_VALUES);
    }
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  VecZeroEntries(b);
  if (!rank) {
    value = 6.66134e-16;
    VecSetValue(b,0,value,INSERT_VALUES);
    value = -value;
    VecSetValue(b,ndof,value,INSERT_VALUES);
  }
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Ac.txt",&viewer);
  MatView(A,viewer);
  PetscViewerDestroy(&viewer);

  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"bc.txt",&viewer);
  VecView(b,viewer);
  PetscViewerDestroy(&viewer);

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPBCGS);
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp,A,A);
  KSPSolve(ksp,b,x);

  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"xc.txt",&viewer);
  VecView(x,viewer);
  PetscViewerDestroy(&viewer);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
  DMDestroy(&da);
  KSPDestroy(&ksp);
  PetscFinalize();

  if (!rank) printf("End of C test program\n");

}
