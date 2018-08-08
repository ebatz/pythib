# Python binding for TwoHadronsInBox

## Introduction

This package contains Python3 bindings for the [TwoHadronsInBox](https://github.com/ebatz/TwoHadronsInBox) library used for the calculation of scattering amplitudes from Lattice QCD spectra.

## Setup
### Prerequisites

The following libraries need to be accessible:

* [TwoHadronsInBox](https://github.com/ebatz/TwoHadronsInBox/tree/qSqDependence) built as a shared library using `make lib`

   **Note**: Currently the particular branch of TwoHadronsInBox linked above has to be used until those changes have been merged back into the master branch. The branch in question facilitates K matrix parametrizations that need access to both the center-of-mass energy and scattering momentum. This functionality is not currently part of the master branch.
* [Pybind11](https://github.com/pybind/pybind11)

### Installation

Only manual installation is available at this point. Adjust the paths to pybind11 as well as TwoHadronsInBox as necessary in build.sh. The build script then produces a Python module that can be imported and used by any Python3 code.

**Note**: The TwoHadronsInBox library is linked dynamically when the Python module is imported. Therefore you need to make sure that `libBox.so` created by TwoHadronsInBox is accessible to the linker either by putting it into one of the standard search paths, or by adding the path to the shared library to LD\_LIBRARY\_PATH.
