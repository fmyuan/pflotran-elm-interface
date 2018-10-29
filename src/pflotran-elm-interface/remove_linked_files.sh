#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    if [ -e "$file" ]; then
        rm -f ${file##*/}
    fi
done

