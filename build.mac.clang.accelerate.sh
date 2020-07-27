#!/usr/bin/env bash

#PYBIND_PATH=/path/to/pybind11/
# NOTE - we assume that the src code for pythib is at the same level in your
#        path as TwoHadronsInBox
TWOHAD_PATH=`pwd | sed "s/pythib/TwoHadronsInBox/"`

mac_fix="-undefined dynamic_lookup"
mac_lapack="-framework Accelerate"

c++ -O3 -Wall -shared -std=c++11 $mac_fix  -fPIC \
    -I${PYBIND_PATH}/include/ `python3 -m pybind11 --includes` \
    -I${TWOHAD_PATH}/source/ \
    -L${TWOHAD_PATH}/build/ BMatrix.cc \
    -o BMat`python3-config --extension-suffix` -lBox $mac_lapack

BMat=`ls -lrt | grep BMat | tail -1 | awk '{print $NF}'`
echo "BMat $BMat"

echo "fixing libpath in $BMat"
echo "install_name_tool -change libBox.so $TWOHAD_PATH/build/libBox.so $BMat"

install_name_tool -change libBox.so $TWOHAD_PATH/build/libBox.so $BMat

echo ""
echo "Check library"
echo "otool -L $BMat"
otool -L $BMat
