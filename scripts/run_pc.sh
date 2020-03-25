#!/bin/bash


BINLIB=$1
PATHMAIN=$2
SCRIPTPATH=$3
class=120

export OMP_NUM_THREADS=6

#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver $PATHMAIN $class 0 > nasoq_fixed_all.csv

#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver $PATHMAIN $class 2 > nasoq_tuned_all.csv

#bash $SCRIPTPATH/run_all_qps_reslib.sh $BINLIB/mosek_demo $PATHMAIN $class  > mosek_all.csv

#bash $SCRIPTPATH/run_all_qps_gurobi.sh $BINLIB/gurobi_demo $PATHMAIN 9  > gurobi_few.csv

#bash $SCRIPTPATH/run_all_qps_lib.sh $BINLIB/osqp_demo $PATHMAIN $class 0 > osqp_nopolished_rest.csv

#bash $SCRIPTPATH/run_all_qps_lib.sh $BINLIB/osqp_demo $PATHMAIN $class 1 > osqp_polished_rest.csv

#bash $SCRIPTPATH/run_all_qps_reslib.sh $BINLIB/ql_demo $PATHMAIN $class  > ql_all_rest.csv


### Update downdate vs. scratch
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_ldlparsy $PATHMAIN 8 0 > nasoq_fixed_scratch_all.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_mkl_sameiter $PATHMAIN 8 0 > nasoq_fixed_MKL_all.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_mkl_default $PATHMAIN 8 0 > nasoq_fixed_MKL_default_all.csv


### Range space vs. full space
#bash $SCRIPTPATH/run_all_qps_reslib.sh $BINLIB/quadprog_demo $PATHMAIN $class  > quadprog_all_rest.csv


### LDL-ParSy
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/ldl_parsy $PATHMAIN $class 0 > LDL_ParSy_Init.csv
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/MKLCholesky $PATHMAIN $class 0 > MKL_Init.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_init $PATHMAIN $class 0 > nasoq_fixed_all_Init.csv


### CHOLMOD row_mod
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_cholmod $PATHMAIN $class 0 > nasoq_cholmod0_p6.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver $PATHMAIN $class 3 1 9 > nasoq_predet0_p9.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver $PATHMAIN $class 3 1 10 > nasoq_predet0_p10.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver $PATHMAIN $class 3 1 8 > nasoq_predet0_p8.csv

### LDL Tuning
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/ldl_parsy_one $PATHMAIN $class 0 > LDL_ParSy_Init_one.csv
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/ldl_parsy_zero $PATHMAIN $class 0 > LDL_ParSy_Init_zero.csv
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/ldl_parsy_min1 $PATHMAIN $class 0 > LDL_ParSy_Init_min1.csv
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/ldl_parsy_min3 $PATHMAIN $class 0 > LDL_ParSy_Init_min3.csv
#bash $SCRIPTPATH/run_all_linsolve.sh $BINLIB/ldl_parsy_min3 $PATHMAIN $class 0 > LDL_ParSy_Init_min3.csv

### LDL Suite
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./linear_solve_wrapper0  ~/UFDB/indef/kkt/ > ldl_suite_l0.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./linear_solve_wrappermin1  ~/UFDB/indef/kkt/ > ldl_suite_lm1.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./linear_solve_wrappermin2  ~/UFDB/indef/kkt/ > ldl_suite_lm2.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./linear_solve_wrappermin3  ~/UFDB/indef/kkt/ > ldl_suite_lm3.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./linear_solve_wrapperpls1  ~/UFDB/indef/kkt/ > ldl_suite_lp1.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./MKLCholesky_metis4  ~/UFDB/indef/kkt/ > mkl_ldl_suite_metis4.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./MKLCholesky_metis5  ~/UFDB/indef/kkt/ > mkl_ldl_suite_metis5.csv
#bash ~/development/pars/scripts/run_lbl_suite.sh  ./linear_solve_wrapper  ~/UFDB/indef/kkt/ > ldl_suite_tuned.csv

### NASOQ uning
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_fixed_min1 $PATHMAIN $class 0 > nasoq_fixed_min1.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_fixed_min2 $PATHMAIN $class 0 > nasoq_fixed_min2.csv
#bash $SCRIPTPATH/run_all_qps.sh $BINLIB/QPSolver_fixed_pls1 $PATHMAIN $class 0 > nasoq_fixed_pls1.csv


