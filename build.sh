#!/usr/bin/env bash

PYBIND_PATH=/path/to/pybind11/
TWOHAD_PATH=/path/to/TwoHadronsInBox/

g++ -O3 -Wall -shared -std=c++11 -fPIC -I${PYBIND_PATH}/include/ `python3-config --includes` -I${TWOHAD_PATH}/source/ -L${TWOHAD_PATH}/build/ BMatrix.cc -o BMat`python3-config --extension-suffix` -lBox -llapack
