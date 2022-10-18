#!/bin/sh

cd src/pflotran
make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 pflotran

