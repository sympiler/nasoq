#!/bin/sh

export MKLROOT=/home/kazem/programs/intel/compilers_and_libraries_2018.2.199/linux/mkl/
export SUITEROOT=/home/kazem/programs/SuiteSparse/
export METISROOT=/home/kazem/programs/metis-5.1.0/build/Linux-x86_64/

mkdir build
cd build
make clean
cmake -DCMAKE_BUILD_TYPE=Release ..
make


## runing a small QP
out_iter=2
in_iter=2
reg_diag=10
tol=15
eps=6
sol_mod=0 # NASOQ-Fixed

BINLIB=./NASOQ
PATHQP=../data/
OUT=../out.csv

rm -f $OUT
# 1 means the type of inputs which expects equality and inequality to be passed separately.
echo "Tool Name,Problem Name,Hessian dim,Hessian NNZ,# of Eq Constraints,Eq Constraint NNZ,# of Ineq Const,Ineq Constraint NNZ,# of Threads,eps_abs,Outer GMRES Iter,Inner GMRES Iter,GMRES Tol,Diagonal Pert,Status,# of Iterations,Time (s),Active-set Size,Constraint Satisfaction Inf,Residual Lagrangian inf,Primal Obj,Dual Obj,Obj Value,Non-negativity Inf,Complementarity Inf," > $OUT

$BINLIB 1 $PATHQP/multiple_contacts_hP.mtx $PATHQP/multiple_contacts_q none none $PATHQP/multiple_contacts_A.mtx $PATHQP/multiple_contacts_l $reg_diag $out_iter $in_iter $eps $tol $sol_mod >> $OUT
echo "">> $OUT
$BINLIB 1 $PATHQP/last_qp_A.txt_A  $PATHQP/last_qp_d.txt  none none $PATHQP/last_qp_C.txt_C  $PATHQP/last_qp_d.txt $reg_diag $out_iter $in_iter $eps $tol>> $OUT
echo "">> $OUT
$BINLIB 1 $PATHQP/osqp_failure_A.txt_A  $PATHQP/osqp_failure_b.txt  none none $PATHQP/osqp_failure_C.txt_C  $PATHQP/osqp_failure_d.txt $reg_diag $out_iter $in_iter $eps $tol>> $OUT
echo "">> $OUT
$BINLIB 1 $PATHQP/osqp_failure_A2.txt_A  $PATHQP/osqp_failure_b2.txt  none none $PATHQP/osqp_failure_C2.txt_C  $PATHQP/osqp_failure_d2.txt $reg_diag $out_iter $in_iter $eps $tol>> $OUT
echo "">> $OUT
$BINLIB 1 $PATHQP/osqp_failure_A3.txt_A  $PATHQP/osqp_failure_b3.txt  none none $PATHQP/osqp_failure_C3.txt_C  $PATHQP/osqp_failure_d3.txt $reg_diag $out_iter $in_iter $eps $tol>> $OUT
echo "">> $OUT
$BINLIB 1 $PATHQP/small_qp_P_A  $PATHQP/small_qp_q  none none $PATHQP/small_qp_A_C  $PATHQP/small_qp_b $reg_diag $out_iter $in_iter $eps $tol>> $OUT
echo "">> $OUT

echo "Check $OUT  for outputs."
