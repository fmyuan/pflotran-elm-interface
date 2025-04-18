#!/bin/sh

. $PFLOTRAN_DIR/.gitlab/skip_on_error.sh

# configure intel oneapi paths
source /opt/intel/oneapi/setvars.sh
export PATH=/opt/intel/oneapi/compiler/latest/linux/bin/intel64:/opt/intel/oneapi/compiler/latest/linux/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

if [[ -n $SRC_DIR ]]; then
  cd $SRC_DIR
else
  cd $PFLOTRAN_DIR/src/pflotran
fi

make clean
make -j4 pflotran_rxn
# prevent building of pflotran if pflotran_rxn is not built correctly
if [ ! -f pflotran_rxn ]; then
  echo "\n----- pflotran_rxn executable not properly compiled -----\n"
  echo 'failed' > $ARTIFACT_DIR/status
else
  make clean
  make -j4 pflotran
fi
