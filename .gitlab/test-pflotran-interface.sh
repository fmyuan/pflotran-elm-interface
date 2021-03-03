#!/bin/sh

export PFLOTRAN_DIR=$PWD
export SRC_DIR=$PFLOTRAN_DIR/src/clm-pflotran
export ARTIFACT_DIR=/tmp/test-pflotran-interface

cd $SRC_DIR

# Run regression tests
make test
REGRESSION_EXIT_CODE=$?

lcov --capture --directory . --output-file pflotran_coverage.info
genhtml pflotran_coverage.info --output-directory coverage

rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR/regression_tests
cp -R regression_tests/* $ARTIFACT_DIR/regression_tests/.
cp pflotran-tests-* $ARTIFACT_DIR/regression_tests/.
cp -R ./coverage $ARTIFACT_DIR

if [ $REGRESSION_EXIT_CODE -ne 0 ]; then
  echo "Regression tests failed" >&2
  exit 1
else
  exit 0
fi

