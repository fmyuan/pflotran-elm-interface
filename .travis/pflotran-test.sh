#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

if [ $CMAKE_BUILD -eq 0 ]; then

  cd src/pflotran;

  # Run unit tests
  make utest
  UNIT_EXIT_CODE=$?
  if [ $UNIT_EXIT_CODE -ne 0 ]; then
    echo "Unit tests failed" >&2
  fi

  # Run regression tests
  cd ../../regression_tests

  make test
  REGRESSION_EXIT_CODE=$?
  if [ $REGRESSION_EXIT_CODE -ne 0 ]; then
    echo "Regression tests failed" >&2
  fi

  if [ $UNIT_EXIT_CODE -eq 0 ] && [ $REGRESSION_EXIT_CODE -eq 0 ]; then
    exit 0
  else
    exit 1
  fi

else

  cd src/pflotran/build;
  ctest --verbose
  REGRESSION_EXIT_CODE=$?
  echo "REGRESSION_EXIT_CODE = " ${REGRESSION_EXIT_CODE}

  if [ $REGRESSION_EXIT_CODE -ne 0 ]; then
    echo "Regression tests failed" >&2
    exit 1
  else
    exit 0
  fi

fi
