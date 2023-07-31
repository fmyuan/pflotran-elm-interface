#!/bin/sh

cd src/pflotran
make clean
make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1 pflotran

