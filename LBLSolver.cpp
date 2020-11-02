//
// Created by kazem on 11/1/20.
//

#include "LBLSolver.h"
namespace nasoq_lib{
 LBLSolver::LBLSolver(int n, int nzmax, int *Ap, int *Ai, double *Ax, double *rhs) {
  // TODO: instantiate ss_
  auto* A = new nasoq::CSC;
  A->nzmax = nzmax;
  A->ncol = A->nrow = n;
  A->packed = 1;

  A->i = Ai;
  A->p = Ap;
  A->x = Ax;

  ss_ = new nasoq::SolverSettings(A, rhs);

  // Use fixed settings for now
  ss_->ldl_variant = 4;
  ss_->req_ref_iter = 0;
  ss_->solver_mode = 0;
  ss_->reg_diag = pow(10,-9);
 }

 LBLSolver::~LBLSolver() {
  delete ss_;
  //TODO deleting A
 }

 int LBLSolver::symbolic_fact() {
  return ss_->symbolic_analysis();
 }

 int LBLSolver::numeric_fact() {
  return ss_->numerical_factorization();
 }

 double* LBLSolver::solve() {
  sol_ = ss_->solve_only();
  return sol_;
 }

 double* LBLSolver::solve(double *in_rhs) {
  memcpy(ss_->rhs, in_rhs, n_ * sizeof(double));
  sol_ = ss_->solve_only();
  return sol_;
 }

 double* LBLSolver::solution(){
  return sol_;
 }

}