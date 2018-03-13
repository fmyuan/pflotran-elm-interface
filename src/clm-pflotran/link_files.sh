#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    ln -s ${file} ${file##*/}
done

#for file in `ls ../pflotran/*.F`; do
#    ln -s ${file} ${file##*/}
#done

#for file in `ls ../pflotran/*.h`; do
#    ln -s ${file} ${file##*/}
#done

