#!/bin/sh

cd $SRC_DIR

REGRESSION_EXIT_CODE=-999

# Run regression tests
TEST_LOG='test.log'
make test 2>&1 | tee $TEST_LOG
if [ $(grep -c "Failed : \|Errors : " "$TEST_LOG") -ne 0 ]; then
  echo "\n----- Regression tests failed -----\n"
  REGRESSION_EXIT_CODE=1
elif [ $(grep -c " All tests passed." "$TEST_LOG") -ne 0 ]; then
  echo "\n----- Regression tests succeeded -----\n"
  REGRESSION_EXIT_CODE=0
else
  echo "\n----- Regression tests produced unrecognized result -----\n"
fi

lcov --capture --directory . --output-file pflotran_coverage.info
genhtml pflotran_coverage.info --output-directory coverage

rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR/regression_tests
cp -R regression_tests/* $ARTIFACT_DIR/regression_tests/.
cp pflotran-tests-* $ARTIFACT_DIR/regression_tests/.
cp -R ./coverage $ARTIFACT_DIR

if [ $REGRESSION_EXIT_CODE -ne 0 ]; then
  echo 'failed' > $ARTIFACT_DIR/status
else
  echo 'success' > $ARTIFACT_DIR/status
fi

