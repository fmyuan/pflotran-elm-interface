#!/bin/sh

export PFLOTRAN_DIR=$PWD
export ARTIFACT_DIR=/tmp/test-pflotran

cd $PFLOTRAN_DIR/src/pflotran

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

cd $PFLOTRAN_DIR/src/pflotran
lcov --capture --directory . --output-file pflotran_coverage.info
genhtml pflotran_coverage.info --output-directory $PFLOTRAN_DIR/coverage

rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR
cp -R $PFLOTRAN_DIR/regression_tests $ARTIFACT_DIR
cp -R $PFLOTRAN_DIR/coverage $ARTIFACT_DIR

if [ $UNIT_EXIT_CODE -eq 0 ] && [ $REGRESSION_EXIT_CODE -eq 0 ]; then
  exit 0
else
  exit 1
fi

