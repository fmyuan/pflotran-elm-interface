#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    if [ -e "$file" ]; then
        rm -f ${file##*/}
    fi
done

if [ -f pflotran_no_provenance.F90 ]; then
    unlink pflotran_no_provenance.F90
fi

if [ -f pflotran_provenance.F90 ]; then
    unlink pflotran_provenance.F90
fi

