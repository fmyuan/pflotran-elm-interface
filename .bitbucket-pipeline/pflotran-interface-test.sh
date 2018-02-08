#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=linux-gnu

# Run regression tests
cd src/clm-pflotran;

make test

cat *.testlog

