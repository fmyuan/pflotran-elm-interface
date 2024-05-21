#!/bin/sh

. $PFLOTRAN_DIR/.gitlab/skip_on_error.sh

cd $SRC_DIR

# initialize to an unitialized value, not 0
UNIT_EXIT_CODE=-999
REGRESSION_EXIT_CODE=-999

# Check to ensure that the pflotran executable exists as we do not want to
# rebuild it inadvertently below as error checking flags will be missing
if [ ! -f pflotran ]; then
  echo 'The PFLOTRAN executable does not exist for testing.'
  rm -Rf $ARTIFACT_DIR
  mkdir -p $ARTIFACT_DIR
  echo 'failed' > $ARTIFACT_DIR/status
  exit 1
fi

# Run unit tests
UTEST_LOG='utest.log'
make gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1 \
  utest 2>&1 | tee $UTEST_LOG
# catch failed tests
if [ $(grep -c " FAILURES!!!\|failed" "$UTEST_LOG") -ne 0 ]; then
  echo "\n----- Unit tests failed -----\n"
  UNIT_EXIT_CODE=1
elif [ $(grep -c " Error 1\|Error: \|undefined reference" "$UTEST_LOG") -ne 0 ]; then
  echo "\n----- Unit test code failed to compile -----\n"
  UNIT_EXIT_CODE=1
elif [ $(grep -c " OK" "$UTEST_LOG") -ne 0 ]; then
  echo "\n----- Unit tests succeeded -----\n"
  UNIT_EXIT_CODE=0
else
  echo "\n----- Unit tests produced unrecognized result -----\n"
fi

# Run regression tests
RTEST_LOG='rtest.log'
# RUNTIME_ERROR_CHECKING toggles on tests for planned errors to ensure that 
# runtime error checking is enabled.
make RUNTIME_ERROR_CHECKING=1 rtest 2>&1 | tee $RTEST_LOG
if [ $(grep -c "Failed : \|Errors : " "$RTEST_LOG") -ne 0 ]; then
  echo "\n----- Regression tests failed -----\n"
  REGRESSION_EXIT_CODE=1
elif [ $(grep -c " All tests passed." "$RTEST_LOG") -ne 0 ]; then
  echo "\n----- Regression tests succeeded -----\n"
  REGRESSION_EXIT_CODE=0
else
  echo "\n----- Regression tests produced unrecognized result -----\n"
fi

cd $SRC_DIR
# revise coverage threshold coloring
echo $'genhtml_hi_limit = 75\n genhtml_med_limit = 25' > ~/.lcovrc
lcov --capture --directory . --output-file pflotran_coverage.info
genhtml pflotran_coverage.info --output-directory coverage

# Ensure that dependencies are updatable. This catches errors is module labels
# in use statements.
cd $SRC_DIR
echo "\n----- Running dependency update to ensure it is functioning properly -----\n"
DEPENDENCY_UPDATE_LOG='dependency_update.log'
python3 ../python/pflotran_dependencies.py 2>&1 | tee $DEPENDENCY_UPDATE_LOG
DEPENDENCY_UPDATE_EXIT_CODE=$?
if [ $DEPENDENCY_UPDATE_EXIT_CODE -ne 0 ]; then
  echo "\n----- Dependency update failing -----\n"
else
  echo "\n----- Dependency update running properly -----\n"
fi

rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR/logs
cp -R $SRC_DIR/$UTEST_LOG $ARTIFACT_DIR/logs
cp -R $SRC_DIR/$RTEST_LOG $ARTIFACT_DIR/logs
cp -R $SRC_DIR/$DEPENDENCY_UPDATE_LOG $ARTIFACT_DIR/logs
cp -R $PFLOTRAN_DIR/regression_tests $ARTIFACT_DIR
cp -R $SRC_DIR/coverage $ARTIFACT_DIR

if [ $UNIT_EXIT_CODE -eq 0 ] && [ $REGRESSION_EXIT_CODE -eq 0 ] && \
   [ $DEPENDENCY_UPDATE_EXIT_CODE -eq 0 ]; then
  echo 'success' > $ARTIFACT_DIR/status
else
  echo 'failed' > $ARTIFACT_DIR/status
fi

