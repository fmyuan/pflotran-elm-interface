#!/bin/sh

for file in `ls ../pflotran/*.F90`; do
    ln -s ${file} ${file##*/}
done

