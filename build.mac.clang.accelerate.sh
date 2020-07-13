#!/usr/bin/env bash

#PYBIND_PATH=/path/to/pybind11/
TWOHAD_PATH=/Users/walkloud/work/research/c51/x_files/code/TwoHadronsInBox

mac_fix="-undefined dynamic_lookup"
mac_lapac="-framework Accelerate"

c++ -O3 -Wall -shared -std=c++11 $mac_fix  -fPIC \
    -I${PYBIND_PATH}/include/ `python3 -m pybind11 --includes` \
    -I${TWOHAD_PATH}/source/ \
    -L${TWOHAD_PATH}/build/ BMatrix.cc \
    -o BMat`python3-config --extension-suffix` -lBox $mac_lapac

echo "fixing libpath in BMat.cpython-36m-darwin.so"
echo "install_name_tool -change libBox.so $TWOHAD_PATH/build/libBox.so BMat.cpython-36m-darwin.so"

install_name_tool -change libBox.so $TWOHAD_PATH/build/libBox.so BMat.cpython-36m-darwin.so
