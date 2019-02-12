#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    if [ -e "$file" ]; then
        rm -f ${file##*/}
    fi
done

#for file in `ls ../pflotran/*.F`; do
#    rm -f ${file##*/}
#done

#for file in `ls ../pflotran/*.h`; do
#    rm -f ${file##*/}
#done
