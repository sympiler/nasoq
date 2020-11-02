//
// Created by kazem on 11/1/20.
//

#ifndef NASOQ_LBLSOLVER_H
#define NASOQ_LBLSOLVER_H

#include <linear_solver_wrapper.h>

namespace nasoq{
class LBLSolver {
 SolverSettings *ss_;
 double *sol_;
 int n_;
public:
 LBLSolver(int n, int *Ap, int *Ai, double *Ax, double *rhs);
 ~LBLSolver();


 int symbolic_fact();
 int numeric_fact();
 double *solve();
 double *solve(double *);
 double *solution();
};

}


#endif //NASOQ_LBLSOLVER_H
