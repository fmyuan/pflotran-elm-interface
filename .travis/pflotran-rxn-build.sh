#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

if [ $CMAKE_BUILD -eq 0 ]; then
  cd src/pflotran;
  if [ $MINIMAL_BUILD -eq 0]; then
    make -j4 codecov=1 pflotran_rxn;
  else
    make -j4 pflotran_rxn;
  fi
fi

