#!/bin/bash

if [[ $(uname) == Darwin ]]; then
  export SDKROOT="${CONDA_BUILD_SYSROOT}"
fi

cd python
make pip
export plumed_default_kernel=$PREFIX/lib/libplumedKernel$SHLIB_EXT
$PYTHON -m pip install . --no-deps -vv

