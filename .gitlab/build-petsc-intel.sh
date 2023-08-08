#!/bin/sh

# clone and build petsc (specific tag)
git clone --depth 1 -b $PETSC_VERSION https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
./configure PETSC_ARCH=petsc-arch \
--with-cc=icc --COPTFLAGS='-g -O0 -diag-disable=10441' \
--with-cxx=icpc --CXXOPTFLAGS='-g -O0 -diag-disable=10441' \
--with-fc=ifort --FOPTFLAGS='-g -O0 -diag-disable=10441 -diag-disable=5462' \
--with-clanguage=c --with-debug=1 --with-shared-libraries=0 --download-hdf5 --download-metis --download-parmetis --download-hypre --download-hdf5-fortran-bindings=yes --download-fblaslapack \
--download-mpich
make
rm -Rf petsc-arch/externalpackages
