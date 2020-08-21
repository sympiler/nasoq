#!/bin/sh


#export MKLROOT=/home/kazem/programs/intel/compilers_and_libraries_2018.2.199/linux/mkl/
#export SUITEROOT=/home/kazem/programs/SuiteSparse/
#export METISROOT=/home/kazem/programs/metis-5.1.0/build/Linux-x86_64/

mkdir build
cd build
make clean
cmake -DMKL_ROOT_PATH=/home/kazem/programs/intel -DMETIS_ROOT_PATH=/home/kazem/programs/metis-5.1.0/build/Linux-x86_64/  -DCMAKE_BUILD_TYPE=Release ..
## You might use the following line for Mac
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/local/Cellar/gcc/9.1.0/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/9.1.0/bin/g++-9 -DMETIS_ROOT_PATH=/Users/kazem/programs/metis-5.1.0/build/Darwin-x86_64 -G "CodeBlocks - Unix Makefiles" ..  -DCMAKE_C_COMPILER=/usr/local/Cellar/gcc/9.1.0/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/9.1.0/bin/g++-9 -DMETIS_ROOT_PATH=/Users/kazem/programs/metis-5.1.0/build/Darwin-x86_64  -DCMAKE_BUILD_TYPE=Release ..
make


## runing a small QP
out_iter=2
in_iter=2
diag_perturb=10
tol=15
eps=6
sol_mod=0 # NASOQ-Fixed

BINLIB=./NASOQ-BIN
PATHQP=../data/
OUT=../data/out.csv
ACC=-3
VARs="fixed tuned predet"
header=1
rm -f $OUT

for VAR in $VARs; do
 for f in $PATHQP/*.yml; do
  if [ $header -eq 1 ]; then
   $BINLIB -i $f -d 1 -e $ACC -v $VAR >>$OUT
   echo "" >>$OUT
   header=0
  else
   $BINLIB -i $f -e $ACC -v $VAR >>$OUT
   echo "" >>$OUT
  fi
 done
done


echo "Check $OUT  for outputs."
