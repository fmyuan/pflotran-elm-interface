#include <stdio.h>

#include "petscmat.h"

int main(int argc,char* argv[]){

  PetscMPIInt rank;
  PetscMPIInt size;

  PetscInt nx;
  PetscInt ny;
  PetscInt n_local;
  PetscInt i, j;
  PetscInt irow, icol;
  PetscInt global_offset;
  PetscInt d_nnz[8];
  PetscInt o_nnz[8];
  PetscReal values[1];

  Mat A;

  PetscCall(PetscInitialize(&argc,&argv,(char *)0,(char *)0));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  if (!rank) printf("Beginning of C test program\n");

  n_local = 8;
  nx = 4;
  ny = 4;
  values[0] = 1.;

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetType(A,MATBAIJ));
  PetscCall(MatSetSizes(A,n_local,n_local,PETSC_DETERMINE,PETSC_DETERMINE));

  if (!rank) printf("Calculating on and off diagonal block fill\n");

  i=0;
  global_offset = n_local*rank;
  if (rank == 0) {
    d_nnz[i] = 3; o_nnz[i++] = 0;
    d_nnz[i] = 4; o_nnz[i++] = 0;
    d_nnz[i] = 4; o_nnz[i++] = 0;
    d_nnz[i] = 3; o_nnz[i++] = 0;
    d_nnz[i] = 3; o_nnz[i++] = 1;
    d_nnz[i] = 4; o_nnz[i++] = 1;
    d_nnz[i] = 4; o_nnz[i++] = 1;
    d_nnz[i] = 3; o_nnz[i++] = 1;
  }
  else {
    d_nnz[i] = 3; o_nnz[i++] = 1;
    d_nnz[i] = 4; o_nnz[i++] = 1;
    d_nnz[i] = 4; o_nnz[i++] = 1;
    d_nnz[i] = 3; o_nnz[i++] = 1;
    d_nnz[i] = 3; o_nnz[i++] = 0;
    d_nnz[i] = 4; o_nnz[i++] = 0;
    d_nnz[i] = 4; o_nnz[i++] = 0;
    d_nnz[i] = 3; o_nnz[i++] = 0;
  }
 
  for (i=0; i<n_local; i++) printf("process %d local row %d global row %d : d_nnz %d o_nnz %d\n",rank,i,i+global_offset,d_nnz[i],o_nnz[i]);

  PetscCall(MatMPIBAIJSetPreallocation(A,1,PETSC_DEFAULT,d_nnz,
                                       PETSC_DEFAULT,o_nnz));

  for (j=ny/2*rank; j<ny/2*(rank+1); j++) {
    for (i=0; i<nx; i++) {
      irow = i + j*nx;
      icol = irow-nx;
      if (j > 0) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,                                                          INSERT_VALUES));
      icol = irow-1;
      if (i > 0) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,
                                               INSERT_VALUES));
      PetscCall(MatSetValuesBlocked(A,1,&irow,1,&irow,values,INSERT_VALUES));
      icol = irow+1;
      if (i < nx-1) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,
                                                  INSERT_VALUES));
      icol = irow+nx;
      if (j < ny-1) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,
                                                  INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  if (!rank) printf("Finished setting up matrix\n");

  PetscCall(MatZeroEntries(A));
  for (j=ny/2*rank; j<ny/2*(rank+1); j++) {
    for (i=0; i<nx; i++) {
      irow = i + j*nx;
      icol = irow-nx;
      if (j > 0) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,                                                          INSERT_VALUES));
      icol = irow-1;
      if (i > 0) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,
                                               INSERT_VALUES));
      PetscCall(MatSetValuesBlocked(A,1,&irow,1,&irow,values,INSERT_VALUES));
      icol = irow+1;
      if (i < nx-1) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,
                                                  INSERT_VALUES));
      icol = irow+nx;
      if (j < ny-1) PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,
                                                  INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  if (!rank) printf("Finished setting matrix values\n");
  PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));

  if (!rank) printf("Begin extending matrix values for well cells\n");
  PetscCall(MatZeroEntries(A));
  PetscCall(MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
  if (rank == 0) {
    irow = 1; icol = 9;
    PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,INSERT_VALUES));
    irow = 1; icol = 13;
    PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,INSERT_VALUES));
    irow = 5; icol = 13;
    PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,INSERT_VALUES));
    irow = 9; icol = 1;
    PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,INSERT_VALUES));
    irow = 13; icol = 1;
    PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,INSERT_VALUES));
    irow = 13; icol = 5;
    PetscCall(MatSetValuesBlocked(A,1,&irow,1,&icol,values,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE));
  if (!rank) printf("Finished extending matrix values for well cells\n");
  PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));

  if (!rank) printf("Finished setting additional matrix value\n");

  if (!rank) printf("End of C test program\n");

  MPI_Finalize();

}
