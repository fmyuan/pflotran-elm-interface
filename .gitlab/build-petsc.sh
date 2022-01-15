#!/bin/sh

git clone https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
git checkout $PETSC_VERSION
./configure PETSC_ARCH=petsc-arch --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --CFLAGS='-g -O0' --CXXFLAGS='-g -O0' --FFLAGS='-g -O0 -Wno-unused-function' --with-clanguage=c --with-debug=1 --with-shared-libraries=0 --download-hdf5 --download-metis --download-parmetis --download-fblaslapack --download-mpich=http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz --download-hypre --download-hdf5-fortran-bindings=yes
make
rm -Rf petsc-arch/externalpackages
