//
// Created by kazem on 2020-05-18.
//

#ifndef NASOQ_LBL_EIGEN_H
#define NASOQ_LBL_EIGEN_H

#include <nasoq/QP/linear_solver_wrapper.h>


 namespace nasoq
 {
  //  Solving Ax=b
  // Inputs:
  //   H  n by n sparse Hessian matrix **lower triangle only** (see
  //     .triangularView<Eigen::Lower>() )
  //   q  n by 1 vector
  // Outputs:
  //   x  n by 1 solution vector
  // Returns nasoq exit flag
  //
  int linear_solve(
    // Pass inputs by copy so we get non-const and casted data
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A,
    Eigen::Matrix<double,Eigen::Dynamic,1> b,
    Eigen::Matrix<double,Eigen::Dynamic,1> & x){
   assert(A.isApprox(A.triangularView<Eigen::Lower>(),0) &&
          "P should be lower triangular");
   assert(A.isCompressed());
   assert(A.rows()==A.cols());
   assert(A.rows()==b.rows());

   CSC *H = new CSC;
   H->nzmax = A.nonZeros();
   H->ncol= H->nrow=A.rows();
   H->p = A.outerIndexPtr();
   H->i = A.innerIndexPtr();
   H->x = A.valuePtr();
   H->stype=-1;
   // TODO: fix this
   H->xtype=CHOLMOD_REAL;
   H->packed=TRUE;
   H->nz = NULL;
   H->sorted = TRUE;
   int reg_diag = -9;

   SolverSettings *lbl = new SolverSettings(H,b.data());
   lbl->ldl_variant = 4;
   lbl->req_ref_iter = 2;
   lbl->solver_mode = 0;
   lbl->reg_diag = pow(10,reg_diag);
   lbl->symbolic_analysis();
   lbl->numerical_factorization();
   double *sol = lbl->solve_only();
   //lbl->compute_norms();
   //std::cout<<"residual: "<<lbl->res_l1<<"\n";

   x = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
     sol,A.rows(),1);

   // Exitflag TODO
   int exitflag = 0;
   delete lbl;
   delete H;
   delete []sol;
   return exitflag;
  }
 }
#endif //NASOQ_LBL_EIGEN_H
