#!/bin/sh

export PFLOTRAN_DIR=$PWD
export SRC_DIR=$PFLOTRAN_DIR/src/pflotran
export ARTIFACT_DIR=/tmp/test-pflotran

cd $SRC_DIR

# initialize to an unitialized value, not 0
UNIT_EXIT_CODE=-999
REGRESSION_EXIT_CODE=-999

# Check to ensure that the pflotran executable exists as we do not want to
# rebuild it inadvertently below as error checking flags will be missing
if [ ! -f pflotran ]; then
  echo 'The PFLOTRAN executable does not exist for testing.'
  echo 'failed' > $ARTIFACT_DIR/status
  exit 1
fi

# Run unit tests
UTEST_LOG='utest.log'
make gnu_code_coverage=1 utest 2>&1 | tee $UTEST_LOG
if [ $(grep -c " FAILURES!!!\|failed" "$UTEST_LOG") -ne 0 ]; then
  echo "\n----- Unit tests failed -----\n" >&2
  UNIT_EXIT_CODE=1
elif [ $(grep -c " OK" "$UTEST_LOG") -ne 0 ]; then
  echo "\n----- Unit tests succeeded -----\n" >&2
  UNIT_EXIT_CODE=0
else
  echo "\n----- Unit tests produced unrecognized result -----\n" >&2
fi

# Run regression tests
RTEST_LOG='rtest.log'
make rtest 2>&1 | tee $RTEST_LOG
if [ $(grep -c "Failed : \|Errors : " "$RTEST_LOG") -ne 0 ]; then
  echo "\n----- Regression tests failed -----\n" >&2
  REGRESSION_EXIT_CODE=1
elif [ $(grep -c " All tests passed." "$RTEST_LOG") -ne 0 ]; then
  echo "\n----- Regression tests succeeded -----\n" >&2
  REGRESSION_EXIT_CODE=0
else
  echo "\n----- Regression tests produced unrecognized result -----\n" >&2
fi

cd $SRC_DIR
lcov --capture --directory . --output-file pflotran_coverage.info
genhtml pflotran_coverage.info --output-directory coverage

rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR
cp -R $PFLOTRAN_DIR/regression_tests $ARTIFACT_DIR
cp -R $SRC_DIR/coverage $ARTIFACT_DIR

if [ $UNIT_EXIT_CODE -eq 0 ] && [ $REGRESSION_EXIT_CODE -eq 0 ]; then
  echo 'success' > $ARTIFACT_DIR/status
else
  echo 'failed' > $ARTIFACT_DIR/status
fi

