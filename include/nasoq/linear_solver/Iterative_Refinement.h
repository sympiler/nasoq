//
// Created by kazem on 11/23/18.
//

#ifndef PROJECT_ITERATIVE_REFINEMENT_H
#define PROJECT_ITERATIVE_REFINEMENT_H



#include "nasoq/linear_solver/solve_phase.h"
#include "nasoq/common/Norm.h"
#include "nasoq/common/SparseUtils.h"
#include "nasoq/matrixVector/spmv_CSC.h"

namespace nasoq {
 int iterative_refinement(int n, size_t *Ap, int *Ai, double *Ax,
                          size_t *Lp, int *Li, double *Lx, int NNZ,
                          size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
                          double *d_val,
                          double epsilon,
                          double *x_h, double *rhs,
                          double *r,
                          int iter_max,
                          double &BWDError) {
  int num_iter = 0;
  double alpha = -1.0;
  double alp[2] = {-1.0, 0};
  double bet[2] = {1.0, 0};
  double normA, normRHS, normX;
  double normR = 0;
  normRHS = norm_dense(1, n, rhs, 0);
  double normxh = norm_dense(1, n, x_h, 0);
  normA = norm_sparse(n, Ap, Ai, Ax, -1, 0);
  //initial error without refinement
  double normAxb = (normA * normxh + normRHS);
  BWDError = normxh / normAxb;
  do {
   //x = A*x_h
   //spmv_csc(n, Ap, Ai, Ax, x_h, x);
   //r = rhs - A*x_h
   for (int i = 0; i < n; ++i) {
    r[i] = rhs[i];
   }
   spmv_csc_sym_one(n, Ap, Ai, Ax, -1, alp, bet, 1, x_h, r);
   normR = norm_dense(1, n, r, 0);

   //solving LDLT * diff = r for diff
   solve_phase_ldl(n, d_val, r, col2sup, sup2col,
                   Lp, Li, Lx, Li_ptr, supNo, NNZ);

   //xh = xh + diff
   add_vec(n, r, 1, x_h);
#if 0
   for (int j = 0; j < n; ++j) {
    if(r[j] > epsilon){
     std::cout<<": "<<j<<";"<<r[j]<<"\n";
    }
   }
   std::cout<<"\n";
#endif
   normX = norm_dense(1, n, x_h, 0);
   //Backward error = |Xc| / (|A||Xinit|+|RHS|)
   BWDError = normX / normAxb;
   num_iter++;
  } while (num_iter < iter_max);
  return num_iter;
 }


 int iterative_refinement_ll(int n, size_t *Ap, int *Ai, double *Ax,
                             size_t *Lp, int *Li, double *Lx, int NNZ,
                             size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
                             double epsilon,
                             double *x_h, double *rhs,
                             double *r,
                             int iter_max,
                             double &BWDError) {
  int num_iter = 0;
  double alpha = -1.0;
  double alp[2] = {-1.0, 0};
  double bet[2] = {1.0, 0};
  double normA, normRHS, normX;
  double normR = 0;
  normRHS = norm_dense(1, n, rhs, 0);
  double normxh = norm_dense(1, n, x_h, 0);
  normA = norm_sparse(n, Ap, Ai, Ax, -1, 0);
  //initial error without refinement
  double normAxb = (normA * normxh + normRHS);
  BWDError = normxh / normAxb;
  //std::cout<<"Initi E: "<< BWDError<<"\n";
  do {
   //x = A*x_h
   //r = rhs - A*x_h
   for (int i = 0; i < n; ++i) {
    r[i] = rhs[i];
   }
   spmv_csc_sym_one(n, Ap, Ai, Ax, -1, alp, bet, 1, x_h, r);
   normR = norm_dense(1, n, r, 0);

   //solving LDLT * diff = r for diff
   solve_phase_ll_blocked(n, r, col2sup, sup2col,
                          Lp, Li, Lx, Li_ptr, supNo, NNZ);

   //xh = xh + diff
   add_vec(n, r, 1, x_h);
#if 0
   for (int j = 0; j < n; ++j) {
    if(r[j] > epsilon){
     std::cout<<": "<<j<<";"<<r[j]<<"\n";
    }
   }
   std::cout<<"\n";
#endif
   normX = norm_dense(1, n, x_h, 0);
   //Backward error = |Xc| / (|A||Xinit|+|RHS|)
   BWDError = normX / normAxb;
   //std::cout<<"==>" <<normR/normRHS<<"\n";
   num_iter++;
  } while (num_iter < iter_max);

  return num_iter;
 }
}
#endif //PROJECT_ITERATIVE_REFINEMENT_H
