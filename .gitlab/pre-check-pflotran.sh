#!/bin/sh

# remove artifact directory if it exists
rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR/logs


cd $SRC_DIR


PRE_CHECK_TEST_WHITESPACE_FLAG=-999
WHITESPACE_LOG="trailing_whitespace.log"
python3 ../python/remove_trailing_whitespace.py > $WHITESPACE_LOG 2>&1
if [ $(grep -c "No trailing whitespaces found" "$WHITESPACE_LOG") -ne 1 ]; then
  cat $WHITESPACE_LOG
  echo "\n----- Trailing whitespace found. These must be removed before code modifications can be accepted. -----\n"
  PRE_CHECK_TEST_WHITESPACE_FLAG=1
else
  PRE_CHECK_TEST_WHITESPACE_FLAG=0
fi

PRE_CHECK_TEST_TAB_FLAG=-999
TAB_LOG="tabs.log"
python3 ../python/remove_tabs.py > $TAB_LOG 2>&1
if [ $(grep -c "No tab characters found" "$TAB_LOG") -ne 1 ]; then
  cat $TAB_LOG
  echo "\n----- Tab found. These must be removed before code modifications can be accepted. -----\n"
  PRE_CHECK_TEST_TAB_FLAG=1
else
  PRE_CHECK_TEST_TAB_FLAG=0
fi

cp $SRC_DIR/$WHITESPACE_LOG $ARTIFACT_DIR/logs
cp $SRC_DIR/$TAB_LOG $ARTIFACT_DIR/logs

if [ $PRE_CHECK_TEST_WHITESPACE_FLAG -eq 0 ] && \
   [ $PRE_CHECK_TEST_TAB_FLAG -eq 0 ]; then
  echo "\n----- Pre-check tests passed. -----\n"
  echo 'success' > $ARTIFACT_DIR/status
else
  echo "\n----- Pre-check tests failed. -----\n"
  echo 'failed' > $ARTIFACT_DIR/status
fi

