#!/bin/sh

cd src/clm-pflotran

# Run regression tests
make test
REGRESSION_EXIT_CODE=$?
if [ $REGRESSION_EXIT_CODE -ne 0 ]; then
  echo "Regression tests failed" >&2
  exit 1
else
  exit 0
fi

