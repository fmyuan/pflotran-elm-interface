#!/bin/sh

cd src/pflotran
make -j4 codecov=1 pflotran

