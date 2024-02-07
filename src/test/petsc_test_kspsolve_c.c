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

  PetscCall(PetscInitialize(&argc,&argv,(char *)0,(char *)0));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  if (!rank) printf("Beginning of C test program\n");

  ncell = 3;
  ndof = 11;

  PetscCall(DMDACreate3d(PETSC_COMM_WORLD,bt,bt,bt,stype,
               ncell,1,1, 
               PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 
               ndof,1,NULL,NULL,NULL,&da));

  PetscCall(DMSetUp(da));
  PetscCall(DMDAGetCorners(da,&offset,&dummy,&dummy,&ncell_local,&dummy,&dummy));
  PetscCall(DMCreateGlobalVector(da,&x));
  PetscCall(VecDuplicate(x,&b));
  PetscCall(DMSetMatType(da,MATBAIJ));
  PetscCall(DMCreateMatrix(da,&A));

  for (i=offset*ndof;i<(offset+ncell_local)*ndof;i++) {
    value = 1333.33;
    PetscCall(MatSetValue(A,i,i,value,INSERT_VALUES));
    if (i > ndof-1) {
      value = -7.5e-7;
      PetscCall(MatSetValue(A,i,i-ndof,value,INSERT_VALUES));
    }
    if (i < (ncell-1)*ndof) {
      value = -7.5e-7;
      PetscCall(MatSetValue(A,i,i+ndof,value,INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  PetscCall(VecZeroEntries(b));
  if (!rank) {
    value = 6.66134e-16;
    PetscCall(VecSetValue(b,0,value,INSERT_VALUES));
    value = -value;
    PetscCall(VecSetValue(b,ndof,value,INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));

  PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Ac.txt",&viewer));
  PetscCall(MatView(A,viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"bc.txt",&viewer));
  PetscCall(VecView(b,viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
  PetscCall(KSPSetType(ksp,KSPBCGS));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetOperators(ksp,A,A));
  PetscCall(KSPSolve(ksp,b,x));

  PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"xc.txt",&viewer));
  PetscCall(VecView(x,viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(DMDestroy(&da));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(PetscFinalize());

  if (!rank) printf("End of C test program\n");

}
