#!/bin/sh

# remove artifact directory if it exists
rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR/logs

PRE_CHECK_TEST_CODE=-999

cd $SRC_DIR

WHITESPACE_LOG="trailing_whitespace.log"
python3 ../python/remove_trailing_whitespace.py > $WHITESPACE_LOG 2>&1
if [ $(grep -c "No trailing whitespaces found" "$WHITESPACE_LOG") -ne 1 ]; then
  cat $WHITESPACE_LOG
  echo "\n----- Trailing whitespace found. These must be removed before code modifications can be accepted. -----\n"
  echo "\n----- Pre-check tests failed. -----\n"
  PRE_CHECK_TEST_CODE=1
else
  echo "\n----- Pre-check tests passed. -----\n"
  PRE_CHECK_TEST_CODE=0
fi

cp $SRC_DIR/$WHITESPACE_LOG $ARTIFACT_DIR/logs

if [ $PRE_CHECK_TEST_CODE -eq 0 ]; then
  echo 'success' > $ARTIFACT_DIR/status
else
  echo 'failed' > $ARTIFACT_DIR/status
fi

