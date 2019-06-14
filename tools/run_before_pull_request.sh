#!/bin/sh
#
# This script builds all versions of PFLOTRAN tested by automated testing
#

echo "\nStart of test for pull request....\n"

cd ..
PFLOTRAN_DIR=`pwd`

GNU_PFLOTRAN_MAKE_SUCCESS=false
GNU_PFLOTRAN_TEST_SUCCESS=false
GNU_PFLOTRAN_RXN_MAKE_SUCCESS=false
#CMAKE_PFLOTRAN_MAKE_SUCCESS=false
#CMAKE_PFLOTRAN_TEST_SUCCESS=false
GNU_INTERFACE_MAKE_SUCCESS=false
GNU_INTERFACE_TEST_SUCCESS=false

# GNU PFLOTRAN ----------------------------------------------------------------

echo "GNU Make PFLOTRAN"

# build PFLOTRAN
MAKE_LOG="make.log"
TEST_LOG="test.log"
cd $PFLOTRAN_DIR/src/pflotran
make clean > /dev/null 2>&1
make pflotran > $MAKE_LOG 2>&1

if [ -e "pflotran" ] ; then
  GNU_PFLOTRAN_MAKE_SUCCESS=true
  echo "  Build passed."

# test PFLOTRAN
  make clean-tests > /dev/null 2>&1
  make test > $TEST_LOG 2>&1

  if [ -e "$TEST_LOG" ] ; then
    # We cannot use Fail or Failure as these are included in the names of
    # several unit test filenames.  The grep catches these when built.
    FAILURES=$(grep -c "Failures" "$TEST_LOG")
    FAILED=$(grep -c "Failed" "$TEST_LOG")
    NUM_FAIL=`expr $FAILURES + $FAILED`
    if [ $NUM_FAIL -gt "0" ] ; then
      echo "  Tests failed.      <--------------------------------PROBLEMS!!!"
    else
      GNU_PFLOTRAN_TEST_SUCCESS=true
      echo "  Tests passed."
    fi
  fi
else
  MYDIR=`pwd`
  echo "  Build failed. See $MYDIR/$MAKE_LOG."
fi

# GNU PFLOTRAN_RXN ------------------------------------------------------------

echo "GNU Make PFLOTRAN_RXN"

# build PFLOTRAN_RXN
MAKE_LOG="pflotran_rxn_make.log"
cd $PFLOTRAN_DIR/src/pflotran
make clean > /dev/null 2>&1
make pflotran_rxn > $MAKE_LOG 2>&1

if [ -e "pflotran_rxn" ] ; then
  GNU_PFLOTRAN_RXN_MAKE_SUCCESS=true
  echo "  Build passed."

# no need to test at this point
else
  MYDIR=`pwd`
  echo "  Build failed. See $MYDIR/$MAKE_LOG."
fi

# GNU PFLOTRAN_INTERFACE ------------------------------------------------------

if [ $GNU_PFLOTRAN_TEST_SUCCESS = true ] ; then

echo "GNU Make PFLOTRAN_INTERFACE"

# build PFLOTRAN_INTERFACE interface
cd $PFLOTRAN_DIR/src/clm-pflotran
./remove_linked_files.sh
./link_files.sh
make clean > /dev/null 2>&1
make pflotran_interface > $MAKE_LOG 2>&1

if [ -e "pflotran_interface" ] ; then
  GNU_INTERFACE_MAKE_SUCCESS=true
  echo "  Build passed."

# test PFLOTRAN_INTERFACE interface
  make clean-tests > /dev/null 2>&1
  make test > $TEST_LOG 2>&1

  if [ -e "$TEST_LOG" ] ; then
    PASSED=$(grep -c "Fail" "$TEST_LOG")
    if [ $PASSED -gt "0" ] ; then
      echo "  Tests failed.      <--------------------------------PROBLEMS!!!"
    else
      echo "  Tests passed."
      GNU_INTERFACE_TEST_SUCCESS=true
    fi
  fi
else
  MYDIR=`pwd`
  echo "  Build failed. See $MYDIR/$MAKE_LOG."
fi

fi # GNU_PFLOTRAN_TEST_SUCCESS

## CMake PFLOTRAN --------------------------------------------------------------
#
#if [ $GNU_INTERFACE_TEST_SUCCESS = true ] ; then
#
#echo "CMake PFLOTRAN"
#
## build PFLOTRAN
#cd $PFLOTRAN_DIR/src/pflotran
#mkdir -p build
#cd build
#rm -Rf *
#cmake ../ -DBUILD_SHARED_LIBS=Off > $MAKE_LOG 2>&1
#make VERBOSE=1 >> $MAKE_LOG 2>&1
#
#if [ -e "pflotran.exe" ] ; then
#  CMAKE_PFLOTRAN_MAKE_SUCCESS=true
#  echo "  Build passed."
#
## test PFLOTRAN
#  ctest --verbose > $TEST_LOG 2>&1
#  if [ -e "$TEST_LOG" ] ; then
#    PASSED=$(grep -c "Fail" "$TEST_LOG")
#    if [ $PASSED -gt "0" ] ; then
#      echo "  Tests failed.      <--------------------------------PROBLEMS!!!"
#    else
#      CMAKE_PFLOTRAN_TEST_SUCCESS=true
#      echo "  Tests passed."
#    fi
#  fi
#else
#  MYDIR=`pwd`
#  echo "  Build failed. See $MYDIR/$MAKE_LOG."
#fi
echo "\nEnd of test for pull request....\n"

#fi # GNU_INTERFACE_TEST_SUCCESS

if [ $GNU_PFLOTRAN_MAKE_SUCCESS = true ] && 
   [ $GNU_PFLOTRAN_TEST_SUCCESS = true ] && 
   [ $GNU_PFLOTRAN_RXN_MAKE_SUCCESS = true ] && 
#   [ $CMAKE_PFLOTRAN_MAKE_SUCCESS = true ] &&
#   [ $CMAKE_PFLOTRAN_TEST_SUCCESS = true ] && 
   [ $GNU_INTERFACE_MAKE_SUCCESS = true ] && 
   [ $GNU_INTERFACE_TEST_SUCCESS = true ] ; then
  echo "You may submit your pull request to bitbucket.org/pflotran/pflotran.\n"
else
  echo "Please fix errors before submitting your pull request. \n"
fi
