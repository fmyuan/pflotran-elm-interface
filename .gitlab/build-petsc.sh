#!/bin/sh

MPICH_DIR=/scratch/mpich-4.1
MPICH_INSTALL_DIR=$MPICH_DIR/install

# clone and build petsc (specific tag)
git clone --depth 1 -b $PETSC_VERSION https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
./configure PETSC_ARCH=petsc-arch \
--with-cc=$MPICH_INSTALL_DIR/bin/mpicc \
--with-cxx=$MPICH_INSTALL_DIR/bin/mpicxx \
--with-fc=$MPICH_INSTALL_DIR/bin/mpif90 \
--COPTFLAGS='-g -O0' --CXXOPTFLAGS='-g -O0' --FOPTFLAGS='-g -O0 -Wno-unused-function' --with-clanguage=c --with-debug=1 --with-shared-libraries=0 --download-hdf5 --download-metis --download-parmetis --download-fblaslapack --download-hypre --download-hdf5-fortran-bindings=yes
make
rm -Rf petsc-arch/externalpackages
