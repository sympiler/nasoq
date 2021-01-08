
This is an example of using NASOQ and LBL out of source tree 
thus it needs nasoq repo to be cloned somewhere that is defined by 
`NASOQ_ROOT` cmake option. If `NASOQ_ROOT` is not defined, the parent 
directory will be used, i.e., `../`. 

You need to first clone NASOQ where you like (`NASOQ_ROOT`)
```
git clone https://github.com/sympiler/nasoq
```

Then you can build project by emitting following commands:
```bash
mkdir build
cd build
cmake -DMKL_ROOT_PATH=path/to/intel -DMETIS_ROOT_PATH=path/to//metis-5.1.0/build/Linux-x86_64/ -DNASOQ_ROOT=path/ended/with/nasoq/  -DCMAKE_BUILD_TYPE=Release ..
make
```

You can also add the path to Suitesparse by adding 
`-DSUITE_ROOT_PATH=path/to/suitesparse`

