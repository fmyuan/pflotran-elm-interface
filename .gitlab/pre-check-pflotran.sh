#!/bin/sh

if [ -n "$ARTIFACT_DIR" ]; then
  LOG_DIR=$ARTIFACT_DIR/logs
  # remove artifact directory if it exists
  rm -Rf $ARTIFACT_DIR
  mkdir -p $LOG_DIR
  cd $SRC_DIR
else
  ARTIFACT_DIR=.
  LOG_DIR=.
  cd $PFLOTRAN_DIR/src/pflotran
fi

printf '  Running tests:\n'
printf '    Checking for trailing whitespace...'
PRE_CHECK_TEST_WHITESPACE_FLAG=-999
WHITESPACE_LOG="$LOG_DIR/trailing_whitespace.log"
python3 ../python/remove_trailing_whitespace.py > $WHITESPACE_LOG 2>&1
if [ $(grep -c "No trailing whitespaces found" "$WHITESPACE_LOG") -ne 1 ]; then
  cat $WHITESPACE_LOG
  echo "\n----- Trailing whitespace found. These must be removed before code modifications can be accepted. -----\n"
  PRE_CHECK_TEST_WHITESPACE_FLAG=1
else
  PRE_CHECK_TEST_WHITESPACE_FLAG=0
  printf 'passed.\n'
fi

printf '    Checking for tabs...'
PRE_CHECK_TEST_TAB_FLAG=-999
TAB_LOG="$LOG_DIR/tabs.log"
python3 ../python/remove_tabs.py > $TAB_LOG 2>&1
if [ $(grep -c "No tab characters found" "$TAB_LOG") -ne 1 ]; then
  cat $TAB_LOG
  echo "\n----- Tab found. These must be removed before code modifications can be accepted. -----\n"
  PRE_CHECK_TEST_TAB_FLAG=1
else
  printf 'passed.\n'
  PRE_CHECK_TEST_TAB_FLAG=0
fi

printf '    Testing dependency generation scripts...'
PRE_CHECK_TEST_DEPENDENCIES_FLAG=-999
DEPENDENCIES_LOG="$LOG_DIR/pflotran_dependencies.log"
python3 ../python/pflotran_dependencies.py test > $DEPENDENCIES_LOG 2>&1
if [ $(grep -c "ERROR:" "$DEPENDENCIES_LOG") -ne 0 ]; then
  cat $DEPENDENCIES_LOG
  echo '\n----- Compilation of PFLOTRAN dependencies failed. This may be due to module names having the correct case in "use" statements. This must be fixed before code modifications can be accepted. -----\n\n'
  PRE_CHECK_TEST_DEPENDENCIES_FLAG=1
else
  printf 'passed.\n'
  PRE_CHECK_TEST_DEPENDENCIES_FLAG=0
fi

printf '    Checking for excessively long source code line lengths...'
PRE_CHECK_TEST_LINE_LENGTH=-999
LINE_LENGTH_LOG="$LOG_DIR/line_length.log"
python3 ../python/check_line_length.py test > $LINE_LENGTH_LOG 2>&1
if [ $(grep -c "No long lines found" "$LINE_LENGTH_LOG") -ne 1 ]; then
  cat $LINE_LENGTH_LOG
  echo '\n----- The source code has lines longer than 132 characters, which violates the PFLOTRAN coding standards. Please reduce to less than 132 characters. This must be fixed before code modifications can be accepted. -----\n'
  PRE_CHECK_TEST_LINE_LENGTH=1
else
  printf 'passed.\n'
  PRE_CHECK_TEST_LINE_LENGTH=0
fi

if [ $PRE_CHECK_TEST_WHITESPACE_FLAG -eq 0 ] && \
   [ $PRE_CHECK_TEST_TAB_FLAG -eq 0 ] && \
   [ $PRE_CHECK_TEST_DEPENDENCIES_FLAG -eq 0 ] && \
   [ $PRE_CHECK_TEST_LINE_LENGTH -eq 0 ]; then
  echo "\n----- Pre-check tests passed. -----\n"
  echo 'success' > $ARTIFACT_DIR/status
else
  echo "\n----- Pre-check tests failed. -----\n"
  echo 'failed' > $ARTIFACT_DIR/status
fi

