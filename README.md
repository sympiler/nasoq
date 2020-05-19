
# NASOQ: Numerically Accurate Sparsity Oriented QP Solver

PLEASE KEEP THE CODE CONFIDENTIAL SINCE IT IS NOT PUBLISHED.

Copyright 2020 Kazem Cheshmi

## Installation
### Library requirements
MKL Pardiso (BLAS) and METIS.

SuiteSparse is an optional dependency which is required when AMD
permutation is enabled.
SuiteSparse includes METIS so if SuiteSparse is installed,
installing METIS is not necessary.

### Building the project
Given that MKL Pardiso and METIS are installed, install NASOQ using
following steps:
```
mkdir build
cd build
cmake -DMKL_ROOT_PATH=path/to/intel -DMETIS_ROOT_PATH=path/to//metis-5.1.0/build/Linux-x86_64/  -DCMAKE_BUILD_TYPE=Release ..
make
```

A quick script for building and running NASOQ is provided in `buildALL.ah`
that you can run it as following:
```bash
bash buildAll.sh
```


### Using NASOQ as a Library
You can use NASOQ as a header file and call it in your application.
You may look at `nasoq_driver.cpp` as an example.

An Eigen interface and driver is provided in`eigen_interface`.

### Testing a QP example
TBD
