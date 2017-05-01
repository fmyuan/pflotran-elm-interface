#!/bin/sh

git clone https://bitbucket.org/petsc/petsc.git

PETSC_GIT_HASH=`cat tools/buildbot/petsc/petsc-git-version.txt`

cd petsc

git checkout ${PETSC_GIT_HASH}

export PETSC_DIR=$PWD
export PETSC_ARCH=petsc-arch

./configure PETSC_ARCH=petsc-arch --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --CFLAGS='-g -O0' --CXXFLAGS='-g -O0' --FFLAGS='-g -O0 -Wno-unused-function' --with-clanguage=c --with-debug=1 --with-shared-libraries=0 --download-hdf5 --download-metis --download-parmetis --download-fblaslapack --download-mpich=http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz

make

