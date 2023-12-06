#include <stdio.h>

#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscis.h"
#include "petscviewer.h"

int main(int argc,char* argv[]){

  PetscMPIInt rank;
  PetscMPIInt size;

  PetscInt ndof;
  PetscInt my_ndof;

  Vec vec;
  PetscViewer viewer;

  PetscCall(PetscInitialize(&argc,&argv,(char *)0,(char *)0));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  if (!rank) printf("Beginning of C test program\n");

  ndof = 100000;
  my_ndof = ((int)(ndof / size)) + (ndof % size) > rank ? 1 : 0;

  PetscCall(VecCreateMPI(PETSC_COMM_WORLD, my_ndof, ndof, &vec));
  PetscCall(VecSet(vec, -999.));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, "vec.bin", FILE_MODE_WRITE, 
                        &viewer));
  PetscCall(PetscViewerBinarySetFlowControl(viewer,2));
  if (!rank) printf("Before VecView\n");
  PetscCall(VecView(vec, viewer));
  if (!rank) printf("After VecView\n");
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(VecDestroy(&vec));

  if (!rank) printf("End of C test program\n");

  MPI_Finalize();

}
