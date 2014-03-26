#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    rm -f ${file##*/}
done

for file in `ls ../pflotran/*.F`; do
    rm -f ${file##*/}
done

for file in `ls ../pflotran/*.h`; do
    rm -f ${file##*/}
done
