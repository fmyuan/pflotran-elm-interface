#!/bin/sh

cd src/clm-pflotran
./link_files.sh
make gnu_code_coverage=1 gnu_runtime_checks=1 pflotran_interface

