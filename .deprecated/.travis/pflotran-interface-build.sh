#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

ls -l ${PETSC_DIR}/${PETSC_ARCH}/lib

cd src/clm-pflotran;
./link_files.sh;
make codedov=1 pflotran_interface;

