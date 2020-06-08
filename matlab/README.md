# NASOQ Matlab wrapper

Simple wrapper for nasoq using libigl and Eigen. There are extra copies, but for
large QPs they should be negligible compared to the nasoq solve.

## Dependencies you must install

MATLAB
cmake
mkl

## Dependencies that cmake will take care of

libigl
eigen

## Installation

From a normal command line terminal:

    mkdir build
    cmake ../ -DCMAKE_BUILD_TYPE=Release
    make

This will create a mex binary in this (`matlab/`) directory.

Then in MATLAB, travel to this folder and issue:

    test_nasoq

## Usage

In MATLAB issue:

    help nasoq
