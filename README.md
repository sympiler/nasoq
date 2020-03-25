
# NASOQ: Numerically Accurate Sparsity Oriented QP Solver

PLEASE KEEP THE CODE CONFIDENTIAL SINCE IT IS NOT PUBLISHED.

Copyright 2020 Kazem Cheshmi

## Installation
### Library requirements
Suitesparse, MKL Pardiso, METIS.

### Building the project

Set the following variables in `buildAll.sh`:
```bash
export MKLROOT <path to MKL>
export SUITEROOT <path to Suitesparse>
export METISROOT <path to METIS> 
```

Then run following script from where you cloned NASOQ:
```bash
bash buildAll.sh
```


### Using NASOQ as a Library
You can use NASOQ as a header file and call it in your application.
You may look at `nasoq_driver.cpp` as an example.

### Testing a QP example
TBD