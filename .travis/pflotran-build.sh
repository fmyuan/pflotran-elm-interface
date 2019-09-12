#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

ls -l ${PETSC_DIR}/${PETSC_ARCH}/lib

if [ $CMAKE_BUILD -eq 0 ]; then
  cd src/pflotran;
  echo 'make -j4 codecov=1 pflotran'
  make -j4 codecov=1 pflotran;
else
  cd src/pflotran
  mkdir build
  cd build
  cmake ../ -DBUILD_SHARED_LIBS=Off
  make VERBOSE=1
fi

