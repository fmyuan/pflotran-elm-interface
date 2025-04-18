#!/bin/sh

. $PFLOTRAN_DIR/.gitlab/skip_on_error.sh

if [[ -n $SRC_DIR ]]; then
  cd $SRC_DIR
else
  cd $PFLOTRAN_DIR/src/pflotran
fi

make clean
make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1 pflotran

