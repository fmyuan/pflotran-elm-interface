#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

ls -l ${PETSC_DIR}/${PETSC_ARCH}/lib

cd src/clm-pflotran;
./link_files.sh;
make pflotran_interface;

