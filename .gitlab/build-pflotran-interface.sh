#!/bin/sh

cd src/clm-pflotran
./link_files.sh
make codecov=1 pflotran_interface

