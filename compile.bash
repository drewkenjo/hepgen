#!/bin/bash

cd $HEPGENDIR/build
#make clean
#cmake ../ -DCMAKE_INSTALL_PREFIX=../install \
#          -DHEPGEN_ENABLE_PYTHON=YES \
#          -DCMAKE_BUILD_TYPE=Release \
#          -DUSE_EXPERIMENTAL=YES \
#          -D${PYTHON_EXECUTABLE}=/apps/python/PRO/bin/python3 \
#          -D${PYTHON_INCLUDE_DIR}=/apps/python/PRO/include/python3.4m/ \
#          -D${PYTHON_LIBRARY}=/apps/python/PRO/lib/libpython3.4m.so
make -j5
make install | grep -v Up-to-date

