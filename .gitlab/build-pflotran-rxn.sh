#!/bin/sh

. $PFLOTRAN_DIR/.gitlab/skip_on_error.sh

cd $SRC_DIR
make clean
make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1 pflotran_rxn
if [ ! -f pflotran_rxn ]; then
  echo "\n----- pflotran_rxn executable not properly compiled -----\n"
  echo 'failed' > $ARTIFACT_DIR/status
fi
