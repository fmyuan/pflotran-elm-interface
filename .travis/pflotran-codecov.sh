#!/bin/sh

if [ $MINIMAL_BUILD -eq 0 ]; then
  bash <(curl -s https://codecov.io/bash)
fi
