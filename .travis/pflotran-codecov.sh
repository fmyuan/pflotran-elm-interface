#!/bin/sh

MINIMAL_BUILD=0
if [ $MINIMAL_BUILD -eq 0 ]; then
  echo 'bash <(curl -s https://codecov.io/bash)'
  bash -c "curl -s https://codecov.io/bash" > curl.stdout
  cat curl.stdout
fi
