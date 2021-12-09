//
// Created by Shujian Qian on 2020-11-17.
//

#ifndef NASOQ_NASOQ_STEP_H
#define NASOQ_NASOQ_STEP_H

#include "nasoq/nasoq.h"

namespace nasoq {
 struct NasoqStep : Nasoq {

  int nxt_active = -1;
  int nxt_drop = -1;
  int to_add = 1;
  int num_violated = 0;
  int compute_norms_per_step = 0;
  double prime_step, dual_step, step;
  solve_type typ_nxt_solve = SOLVE;

  int step_solve_finished;  // 1: finished, 0: not finished
  int step_tuned_iter;

  NasoqStep(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in,
            size_t B_row, size_t B_col, int *Bp, int *Bi,
            double *Bx, double *b_ineq);

  NasoqStep(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in,
            size_t A_size1, size_t A_size2, int *Ap, int *Ai, double *Ax,
            double *a_eq,
            size_t B_size1, size_t B_size2, int *Bp, int *Bi, double *Bx,
            double *b_ineq);

  /**
   * Initialize the QP solver, does Symbolic_QP and initialize_x
   * @return status of the initialization 0:success 1:fail
   */
  int solve_init();

  /**
   * Advance the iterative solve by one step
   * @return nasoq_status of the solve
   */
  int solve_step();

  /**
   * Advance the iterative solve by several steps. Stops early when the solver finishes.
   * Unchecked error: called on a solve that is already finished;
   * @param num_step Number of steps to advance
   * @return nasoq_status of the solve
   */
  int solve_steps(int num_step);

  nasoq_status check_step_status(int status);
 };

}

#endif //NASOQ_NASOQ_STEP_H
