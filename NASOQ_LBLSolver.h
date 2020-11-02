//
// Created by kazem on 11/1/20.
//

#ifndef NASOQ_LBLSOLVER_H
#define NASOQ_LBLSOLVER_H

#include <linear_solver_wrapper.h>

namespace nasoq{
class NASOQ_LBLSolver {
 SolverSettings *ss_;
 double *sol_;
 int n_;
 int nzmax_;
public:
 /// Takes a lower triangular CSC indefinite matrix and a right hand side vector
 /// with size of n
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param rhs
 NASOQ_LBLSolver(int n, int nzmax, int *Ap, int *Ai, double *Ax, double *rhs);
 ~NASOQ_LBLSolver();

 /// Symbolic factorization, happens once and can be reused as long as sparsity does not change
 /// \return
 int symbolic_fact();

 /// Actual factorization happens here
 /// \return
 int numeric_fact();

 /// Solves the system for the RHS given to the constructor.
 /// \return
 double *solve();

 /// Replaces the initial RHS with the given input
 /// \return
 double *solve(double *);
 double *solution();
};

}


#endif //NASOQ_LBLSOLVER_H
