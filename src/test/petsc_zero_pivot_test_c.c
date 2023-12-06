#include <stdio.h>

#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscis.h"
#include "petscviewer.h"

int main(int argc,char* argv[]){

  PetscMPIInt rank;
  PetscMPIInt size;

  Vec x;
  Vec b;
  Mat A;
  KSP ksp;
  PC pc;
  PetscReal value;
  PetscReal tolerance;
  PetscInt i;
  PetscInt n;
  PetscInt offset;
  PetscScalar *x_ptr;

  PetscCall(PetscInitialize(&argc,&argv,(char *)0,(char *)0));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  if (!rank) printf("Beginning of C test program\n");

  n = 2/size;
  PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,n,n,
               PETSC_DETERMINE,PETSC_DETERMINE,
               1,NULL,
               0,NULL,&A));
  PetscCall(VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DETERMINE,&x));
  PetscCall(VecDuplicate(x,&b));
  offset = rank * n;
  value = 1.e-16;
  for(i=0;i<n;i++) {
    PetscCall(MatSetValue(A,i+offset,i+offset,value,INSERT_VALUES));
    PetscCall(VecSetValue(b,i+offset,value,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));

  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
  PetscCall(KSPSetErrorIfNotConverged(ksp,PETSC_TRUE));
  PetscCall(KSPSetOperators(ksp,A,A));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(PCSetFromOptions(pc));
//  KSPSetUp(ksp));

  tolerance = 1.e-20;
  PetscCall(PCFactorSetZeroPivot(pc,tolerance));

//  KSPSetUp(ksp));
  PetscCall(KSPSolve(ksp,b,x));

  PetscCall(VecGetArray(x,&x_ptr));
  if (!rank) printf("These values should be ~1:");
  for(i=0;i<n;i++)
    printf(" %f",x_ptr[i]);
  if (!rank) printf("\n");
  PetscCall(VecRestoreArray(x,&x_ptr));

  PetscCall(KSPDestroy(&ksp));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));

  if (!rank) printf("End of C test program\n");

  MPI_Finalize();

}
