#!/bin/sh

# first build mpich
tar -xzvf mpich-4.1.tar.gz
MPICH_DIR=/scratch/mpich-4.1
MPICH_INSTALL_DIR=$MPICH_DIR/install
cd $MPICH_DIR
# from petsc 3.18 --download-mpich=yes configure  for mpich 4.0.2
./configure --prefix=$MPICH_INSTALL_DIR MAKE=/usr/bin/gmake --libdir=$MPICH_INSTALL_DIR/lib CC=gcc CFLAGS="-g -O0 -fPIC" AR=/usr/bin/ar ARFLAGS=cr CXX=g++ CXXFLAGS="-g -O0 -std=gnu++17" FFLAGS="-g -O0 -Wno-unused-function -fallow-argument-mismatch" FC=gfortran F77=gfortran FCFLAGS="-g -O0 -Wno-unused-function -fallow-argument-mismatch" --disable-shared --with-pm=hydra --disable-java --with-hwloc=embedded --enable-fast=no --enable-error-messages=all --with-device=ch3:sock --enable-g=meminit
make all; make install

# then clone and build petsc
git clone https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
git checkout $PETSC_VERSION
./configure PETSC_ARCH=petsc-arch \
--with-cc=$MPICH_INSTALL_DIR/bin/mpicc \
--with-cxx=$MPICH_INSTALL_DIR/bin/mpicxx \
--with-fc=$MPICH_INSTALL_DIR/bin/mpif90 \
--COPTFLAGS='-g -O0' --CXXOPTFLAGS='-g -O0' --FOPTFLAGS='-g -O0 -Wno-unused-function' --with-clanguage=c --with-debug=1 --with-shared-libraries=0 --download-hdf5 --download-metis --download-parmetis --download-fblaslapack --download-hypre --download-hdf5-fortran-bindings=yes
make
rm -Rf petsc-arch/externalpackages
