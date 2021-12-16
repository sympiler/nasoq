//
// Created by Shujian Qian on 2020-11-17.
//


#include "nasoq/nasoq_step.h"
#include <cstdio>
#include <iostream>

namespace nasoq {

 NasoqStep::NasoqStep(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in,
                          size_t B_row, size_t B_col, int *Bp, int *Bi,
                          double *Bx, double *b_ineq) : Nasoq(H_size, Hp, Hi, Hx, q_in,
                                                              B_row, B_col, Bp, Bi,
                                                              Bx, b_ineq) {}

 NasoqStep::NasoqStep(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in,
                          size_t A_size1, size_t A_size2, int *Ap, int *Ai, double *Ax,
                          double *a_eq,
                          size_t B_size1, size_t B_size2, int *Bp, int *Bi, double *Bx,
                          double *b_ineq) : Nasoq(H_size, Hp, Hi, Hx, q_in,
   A_size1, A_size2, Ap, Ai, Ax,
   a_eq,
   B_size1, B_size2, Bp, Bi, Bx,
   b_ineq) {}

 int NasoqStep::solve_init() {
 ret_val = nasoq_status::NotFinished;
 int status = 1;
 warm_start = 0;
 //print_csc("BT fisrt:\n",BT->ncol,BT->p,BT->i,BT->x);
 qi->tot_st = qi->tic();
 if (scaling > 0) {
  scale_data();
  //print_csc("H scaled\n",H->ncol,H->p,H->i,H->x);
 }
 if (scaling < 0) {
  scale_data_2();
 }
 if (variant == nasoq_mode::Tuned) {//TODO: reuse symbolic QP
  symbolic_QP();
  // the tuned solver loop is now controlled by step_tuned_iter
  initialize_x();

 } else { // Fixed or Auto
  if (variant == nasoq_mode::AUTO) {
   auto_reg_en = 2;
  } else if (variant == nasoq_mode::Fixed) {
   auto_reg_en = 1;
  } else {
   auto_reg_en = 0;
  }
#ifdef CHOLROWMOD
  status = symbolic_QP_cholmod();
    initialize_x_cholmod();
#else
  symbolic_QP();
  status = initialize_x();
#endif
 }

 nxt_active = -1;
 nxt_drop = -1;
 to_add = 1;
 num_violated = 0;
 typ_nxt_solve = SOLVE;
 step_solve_finished = 0;
 step_tuned_iter = 0;

 if (warm_start)
  typ_nxt_solve = REFACTOR;
 return status;
}

 int NasoqStep::solve_step() {
 if (num_iter > max_iter_nas) {
  is_converged = 0;
  step_solve_finished = 1;
  return check_step_status(is_converged);
 }
 if (to_add == 1) {
  num_violated = primal_feasibility(nxt_active);
  if (num_violated == 0) {
   is_converged = 1;
   step_solve_finished = 1;
   return check_step_status(is_converged);
  } else if (nxt_active == -1) {//unbounded!
   std::cout << "All constraints are used\n";
   is_converged = 0;
   step_solve_finished = 1;
   return check_step_status(is_converged);
  }
  update_rhs(nxt_active);
 }
#ifdef CHOLROWMOD
 solve_kkt_cholmod(typ_nxt_solve);
#else
 solve_kkt(typ_nxt_solve);
#endif
 prime_step = primal_step_length(nxt_active);
 dual_step = dual_step_length(nxt_drop);
 step = prime_step > dual_step ? dual_step : prime_step;

 if (step == std::numeric_limits<double>::max()) {
  is_converged = 0;
  step_solve_finished = 1;
  return check_step_status(is_converged);
 } else if (prime_step == std::numeric_limits<double>::max()
            && step == dual_step) { // Dual step
  update_dual(step);
  dual_vars[nxt_active] = dual_vars[nxt_active] + step;//==-step*-1
  to_add = 0;
#ifdef CHOLROWMOD
  edit_kkt_cholmod(to_add, nxt_drop);
#else
  edit_kkt(to_add, nxt_drop);
#endif

  typ_nxt_solve = UPDATE;
  int is_up = update_active_set(to_add, nxt_drop);
  assert(is_up == 1);

 } else { // Full step
  update_primal(step);// solution = step*descent + solution
  update_dual(step);// dual_vars = step*descent + dual_vars
  dual_vars[nxt_active] = dual_vars[nxt_active] + step;//==-step*-1
  if (step == prime_step) { //update
   to_add = 1;
#ifdef CHOLROWMOD
   edit_kkt_cholmod(to_add, nxt_active);
#else
   edit_kkt(to_add, nxt_active);
#endif

   typ_nxt_solve = UPDATE;
   update_active_set(to_add, nxt_active);
  } else if (step == dual_step) { //downdate
   to_add = 0;
#ifdef CHOLROWMOD
   edit_kkt_cholmod(to_add, nxt_drop);
#else
   edit_kkt(to_add, nxt_drop);
#endif
   typ_nxt_solve = UPDATE;
   int is_up = update_active_set(to_add, nxt_drop);
   assert(is_up == 1);
  }
  assert(n_active == active_set.size());
 }
 is_converged = 0;
 return check_step_status(is_converged);
}

 nasoq_status NasoqStep::check_step_status(int status) {

 if (!step_solve_finished) {
  ++num_iter;
  ret_val = nasoq_status::NotFinished;
  return ret_val;
 }

 detect_solver_name();

 if (variant == nasoq_mode::Tuned) {
  if (step_tuned_iter < 4) {
   ++step_tuned_iter;
   initialize_x();
   step_solve_finished = 0;
   ret_val = nasoq_status::NotFinished;
   return ret_val;
  }
 }

 if (scaling > 0) {
  unscale_data();
  unscale_solution();
 }
 if (scaling < 0) {
  unscale_data_2();
  unscale_solution_2();
 }
 qi->tot_end = qi->toc();
 qi->tot = qi->elapsed_time(qi->tot_st, qi->tot_end);

//  ret_val = check_solve_status(status);
 detect_solver_name();
//  return ret_val;

 if(step_solve_finished && status == 0)
  return nasoq_status::Infeasible;
 if(compute_norms_per_step || step_solve_finished){
  compute_objective();
  constraint_sat_norm();
  lagrangian_residual_norm();
  complementarity_norm();
  non_negativity_norm();
 }
 if(step_solve_finished){
  if (lag_res <= eps_abs && cons_sat_norm <= eps_abs &&
      non_negativity_infn <= eps_abs) {//&& complementarity_infn <= eps_abs
   ret_val = nasoq_status::Optimal; //converged
  } else if (cons_sat_norm <= eps_abs) {//low accuracy
   ret_val = nasoq_status::Inaccurate;
  } else {
   ret_val = nasoq_status::NotConverged; //not converged or maybe very inaccurate
  }
 } else{
  ret_val = nasoq_status::NotFinished;
 }
 return ret_val;
}

 int NasoqStep::solve_steps(int num_step) {
 int status = nasoq_status::NotFinished;
 for (int i = 0; i < num_step; i++) {
  status = solve_step();
  if (status != nasoq_status::NotFinished) break;
 }
 return status;
}
}

