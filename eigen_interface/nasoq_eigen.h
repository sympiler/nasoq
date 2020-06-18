//
// Created by kazem on 2020-05-08.
//

#ifndef NASOQ_NASOQ_EIGEN_H
#define NASOQ_NASOQ_EIGEN_H
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <nasoq.h>



 namespace nasoq
 {
  //  minimize        0.5 x' H x + q' x
  // subject to       A x = b; Cx <= d;
  //
  // Inputs:
  //   H  n by n sparse Hessian matrix **lower triangle only** (see
  //     .triangularView<Eigen::Lower>() )
  //   q  n by 1 linear coefficients
  //   A  m1 by n linear equality constraint coefficients
  //   a  m1 by 1 linear equality upper bounds
  //   C  m2 by n linear inequality constraint coefficients
  //   d  m2 by 1 linear inequality upper bounds

  // Outputs:
  //   x  n by 1 solution vector of primal vars
  //   y  m1 by 1 solution vector of primal vars
  //   z  m2 by 1 solution vector of primal vars
  // Returns nasoq exit flag
  //
  int quadprog(
    // Pass inputs by copy so we get non-const and casted data
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> H,
    Eigen::Matrix<double,Eigen::Dynamic,1> q,
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A,
    Eigen::Matrix<double,Eigen::Dynamic,1> b,
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> C,
    Eigen::Matrix<double,Eigen::Dynamic,1> d,
    Eigen::Matrix<double,Eigen::Dynamic,1> & x,
    Eigen::Matrix<double,Eigen::Dynamic,1> & y,
    Eigen::Matrix<double,Eigen::Dynamic,1> & z,
    QPSettings *qs ){
   assert(H.isApprox(H.triangularView<Eigen::Lower>(),0) &&
          "P should be lower triangular");
   assert(H.isCompressed());
   assert(A.isCompressed());
   assert(H.rows()==H.cols());
   assert(H.rows()==q.rows());
   assert(H.cols()==A.cols());
   assert(A.rows()==b.rows());
   assert(C.rows()==d.rows());
   // Exitflag
   int exitflag = 0;
   /// Populate data
   auto nasoq = new Nasoq(H.rows(),H.outerIndexPtr(),H.innerIndexPtr(),H.valuePtr(),q.data(),
                      A.rows(),A.cols(),A.outerIndexPtr(),A.innerIndexPtr(),
                      A.valuePtr(),b.data(),
                      C.rows(),C.cols(),C.outerIndexPtr(),C.innerIndexPtr(),
                      C.valuePtr(),d.data());
   /// Define solver settings if provided
   if(qs){
    nasoq->reg_diag=qs->reg_diag;
    nasoq->eps_rel = qs->eps;
    nasoq->eps_abs = qs->eps;
    nasoq->zero_thresh = qs->zero_thresh;
    nasoq->inner_iter_ref = qs->inner_iter_ref;
    nasoq->outer_iter_ref = qs->inner_iter_ref;
    nasoq->tol_ref = qs->tol_ref;
    nasoq->max_iter = qs->max_iter;
    if (qs->nasoq_mode == "fixed")
     nasoq->nas_mode = Fixed;
    else if (qs->nasoq_mode == "tuned")
     nasoq->nas_mode = Tuned;
    else if(qs->nasoq_mode == "auto")
     nasoq->nas_mode = AUTO;
    else
     nasoq->nas_mode= PREDET;
   }

   /// Solve Problem
   exitflag = nasoq->solve();

   /// Copy output
   x = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
     nasoq->primal_vars,H.rows(),1);
   y = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
     nasoq->dual_vars_eq,A.rows(),1);
   z = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
     nasoq->dual_vars,C.rows(),1);

   /// free memory
   delete nasoq;
   return exitflag;
  }
 }


#endif //NASOQ_NASOQ_EIGEN_H
