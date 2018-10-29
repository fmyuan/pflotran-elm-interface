#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    ln -sf ${file} ${file##*/}
done

if [ -f pflotran.F90 ]; then
    unlink pflotran.F90
fi

if [ -f pflotran_rxn.F90 ]; then
    unlink pflotran_rxn.F90
fi

if [ ! -f pflotran_provenance.F90 ]; then
    mv pflotran_no_provenance.F90 pflotran_provenance.F90
else
    unlink pflotran_no_provenance.F90
fi