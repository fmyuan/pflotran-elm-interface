#!/bin/sh

. $PFLOTRAN_DIR/.gitlab/skip_on_error.sh

if [ -n "$SRC_DIR" ]; then
  echo 'Using SRC_DIR.'
  cd $SRC_DIR
else
  echo 'Using PFLOTRAN_DIR/src/pflotran.'
  cd $PFLOTRAN_DIR/src/pflotran
fi

make clean
make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1 pflotran_rxn
if [ ! -f pflotran_rxn ]; then
  echo "\n----- pflotran_rxn executable not properly compiled -----\n"
  echo 'failed' > $ARTIFACT_DIR/status
fi
