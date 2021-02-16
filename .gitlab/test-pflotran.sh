#!/bin/sh

export PFLOTRAN_DIR=$PWD
export SRC_DIR=$PFLOTRAN_DIR/src/pflotran
export ARTIFACT_DIR=/tmp/test-pflotran

cd $SRC_DIR

# Run unit tests
make codecov=1 utest

UNIT_EXIT_CODE=$?
if [ $UNIT_EXIT_CODE -ne 0 ]; then
  echo "Unit tests failed" >&2
fi

# Run regression tests
cd $PFLOTRAN_DIR/regression_tests

make test
REGRESSION_EXIT_CODE=$?
if [ $REGRESSION_EXIT_CODE -ne 0 ]; then
  echo "Regression tests failed" >&2
fi

cd $SRC_DIR
lcov --capture --directory . --output-file pflotran_coverage.info
genhtml pflotran_coverage.info --output-directory coverage

rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR
cp -R $PFLOTRAN_DIR/regression_tests $ARTIFACT_DIR
cp -R $SRC_DIR/coverage $ARTIFACT_DIR

if [ $UNIT_EXIT_CODE -eq 0 ] && [ $REGRESSION_EXIT_CODE -eq 0 ]; then
  exit 0
else
  exit 1
fi

