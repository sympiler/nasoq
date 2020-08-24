
This is an example of using NASOQ and LBL out of source tree 
thus it needs nasoq repo to be cloned within this
subdirectory, i.e., `nasoq/examples/`

You need to first clone NASOQ where this project is
```
git clone https://github.com/sympiler/nasoq
```

Then you can build project by emitting following commands:
```bash
mkdir build
cd build
cmake -DMKL_ROOT_PATH=path/to/intel -DMETIS_ROOT_PATH=path/to//metis-5.1.0/build/Linux-x86_64/  -DCMAKE_BUILD_TYPE=Release ..
make
```

You can also add the path to Suitesparse by adding 
`-DSUITE_ROOT_PATH=path/to/suitesparse`

