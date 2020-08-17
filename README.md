
# NASOQ: Numerically Accurate Sparsity Oriented QP Solver

NASOQ is a scalable and efficient Quadratic Programming solver that 
obtains solutions for requested accuracies. Visit our website for more details: 
[NASOQ Website](https://nasoq.github.io/)




## Installation
### Library requirements
MKL Pardiso (BLAS), OpenMP and METIS.


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
Upon successful build you should be able to see `data/out.csv` and 
it should be similar to `data/out_correct.csv`.  

For installing on MAc you might need to use GCC so you need to also set the CMAKE compiler flag.

More details are provided in: https://nasoq.github.io/docs/getting-started-nasoq/

### Using NASOQ as a Library
You can use NASOQ as a header file and call it in your application.
More details: https://nasoq.github.io/docs/getting-started-nasoq/


### Testing a QP example
To test a QP example you may also use `NASOQ-BIN` which is a command line interfce for NASOQ.
Some small QP problems are available in `data` folder.
For evaluating NASOQ versus other solvers a separate repository is also provided in:
https://github.com/sympiler/nasoq-benchmarks
More details: https://nasoq.github.io/docs/repository/


Copyright 2020 Kazem Cheshmi
