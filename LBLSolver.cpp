//
// Created by kazem on 11/1/20.
//

#include "LBLSolver.h"
namespace nasoq{
 LBLSolver::LBLSolver(int n, int *Ap, int *Ai, double *Ax, double *rhs) {
  // TODO: instantiate ss_
 }

 LBLSolver::~LBLSolver() {
  delete ss_;
 }

 int LBLSolver::symbolic_fact() {
  return ss_->symbolic_analysis();
 }

 int LBLSolver::numeric_fact() {
  return ss_->numerical_factorization();
 }

 double* LBLSolver::solve() {
  return ss_->solve_only();
 }

 double* LBLSolver::solve(double *in_rhs) {
  memcpy(ss_->rhs, in_rhs, n_ * sizeof(double));
  return ss_->solve_only();
 }

}