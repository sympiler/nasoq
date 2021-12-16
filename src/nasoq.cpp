//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/nasoq.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "nasoq/common/Norm.h"
#include "nasoq/common/transpose_unsym.h"
#include "nasoq/common/Util.h"
#include "nasoq/matrixVector/spmv_CSC.h"
#include "nasoq/QP/nasoq_utils.h"
#include "nasoq/QP/osqp_utils.h"
#include "nasoq/QP/updown_test.h"
#include "nasoq/common/Sym_BLAS.h"

namespace nasoq {
 QPSettings::QPSettings() {
  eps = eps_rel = pow(10, -3);
  eps_primal = eps_dual = eps_slack = eps_nn = eps;
  diag_perturb = zero_thresh = pow(10, -9);
  max_iter_nas = 4000;
  batch_size = 1;
  scaling = 0;
  inner_iter_ref = 0;
  outer_iter_ref = 0; max_iter=0;
  stop_tol = 1e-15;
  nasoq_variant = "Fixed";
 }

 nasoq_config::nasoq_config(int a, int b, double c, double d)
   : inner_iter(a), outer_iter(b), pert_diag(c), stop_tol(d) {}

 Nasoq::Nasoq(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in, size_t B_row, size_t B_col, int *Bp, int *Bi,
              double *Bx, double *b_ineq) {
  H = new CSC;
  H->nzmax = Hp[H_size];
  H->ncol = H->nrow = H_size;
  H->stype = -1;
  H->xtype = CHOLMOD_REAL;
  H->packed = TRUE;
  H->p = Hp;
  H->i = Hi;
  H->x = Hx;
  H->nz = NULL;
  H->sorted = TRUE;
  hessian_size = H->nrow;
  //HT = new CSC;
  q = q_in;

  A = new CSC;
  A->nzmax = 0;
  A->ncol = 0;
  A->nrow = 0;
  A->stype = 1;
  A->xtype = CHOLMOD_REAL;
  A->packed = TRUE;
  A->p = NULL;
  A->i = NULL;
  A->x = NULL;
  A->nz = NULL;
  A->sorted = TRUE;
  AT = new CSC;
  a = NULL;
  eq_const_size = A->nrow;

  B = new CSC;
  B->ncol = B_col;
  B->nrow = B_row;
  B->stype = 0;
  B->xtype = CHOLMOD_REAL;
  B->packed = TRUE;
  B->p = Bp;
  B->i = Bi;
  B->x = Bx;
  B->nz = NULL;
  B->sorted = TRUE;
  BT = new CSC;
  b = b_ineq;
  ineq_const_size = B->nrow;
  if (B->nrow > 0) {
   B->nzmax = Bp[B_col];
  }else{
   B->nzmax=0;
  }

  //sKKT = new CSC;
  //sKKTt = new CSC;
  workspace = new double[B->nrow + 4 * H->ncol];
  primal_vars = new double[H->ncol];
  dual_vars = new double[B->nrow]();
  used_const = new int[B->nrow]();
  rec_length = new double[B->nrow]();
  num_active = 0;
  num_iter = 0;
  warm_start = 0;
  n_active = 0;
  qi = new qp_info;
  default_setting();

#ifdef CHOLROWMOD
  L = NULL;
   kkt = NULL;
   kkt_updated = new CSC;
   cm = &Common;
   cholmod_l_start(cm);
   //cm->final_asis = 1;//TRUE
  // cm->nmethods = 1;
   cm->method[0].ordering = CHOLMOD_METIS;
   k_size = 0;
#endif
 }

 Nasoq::Nasoq(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in, size_t A_size1, size_t A_size2, int *Ap,
              int *Ai, double *Ax, double *a_eq, size_t B_size1, size_t B_size2, int *Bp, int *Bi, double *Bx,
              double *b_ineq) :
   Nasoq(H_size, Hp, Hi, Hx, q_in, B_size1, B_size2, Bp, Bi, Bx, b_ineq) {

  A->ncol = A_size2;
  A->nrow = A_size1;
  A->stype = 1;
  A->xtype = CHOLMOD_REAL;
  A->packed = TRUE;
  A->p = Ap;
  A->i = Ai;
  A->x = Ax;
  A->nz = NULL;
  A->sorted = TRUE;
  a = a_eq;
  if (A->nrow > 0) {
   A->nzmax = Ap[H_size];
   dual_vars_eq = new double[A->nrow]();
  }
  eq_const_size = A->ncol;
 }

 Nasoq::~Nasoq() {
  //TODO
  delete[]workspace;
  //delete []sKKTrhs;
  //delete []kkt_solution;
  if (scaling > 0) {
   delete[]D;
   delete[]E;
   delete[]F;
   delete[]Finv;
   delete[]Einv;
   delete[]Dinv;
  }
  delete[]primal_vars;
  if (A->nrow > 0)
   delete[]dual_vars_eq;
  delete[]dual_vars;
  delete[]used_const;
#ifndef CHOLROWMOD
  delete ss;
#endif
  delete qi;
  if (B->nrow > 0)
   allocateAC(BT, 0, 0, 0, FALSE);
  else
   delete BT;
  if (A->nrow > 0)
   allocateAC(AT, 0, 0, 0, FALSE);
  else
   delete AT;
  //allocateAC(sKKT,0,0,0,FALSE);
  //allocateAC(sKKTt,0,0,0,FALSE);
  delete[]rec_length;
  delete H;
  delete A;
  delete B;
#ifdef CHOLROWMOD
  cholmod_l_free_factor(&L, cm);
   cholmod_l_free_sparse(&kkt, cm);
   //allocateAC(kkt_updated,0,0,0,FALSE);
#endif
 }

 void Nasoq::default_setting() {
  max_iter_nas = 4000;
  diag_perturb = 1e-9;
  zero_thresh = diag_perturb;
  eps_abs = 1e-3;
  eps_rel = eps_abs; // no use
  warm_start = 0;
  batch_size = 1;
  scaling = 0;
  inner_iter_ref = 0;
  outer_iter_ref = 0;
  max_iter=0;
  stop_tol = 1e-15;
  auto_reg_en = 0;
  variant = nasoq_mode::Fixed;
  sol_name = "NASOQ-Fixed";
 }

 int Nasoq::symbolic_QP() {
  int status = 0;
  transpose_unsym(B->nrow, B->ncol, B->p, B->i, B->x,
                  BT->nrow, BT->ncol, BT->p, BT->i, BT->x);
  BT->nzmax = B->nzmax;
  //print_csc("BT:\n",BT->ncol,BT->p,BT->i,BT->x);
  //build_super_kkt();
  //Creates a new solver instance
  double *sKKTrhs_init = new double[hessian_size];
  if (A->nrow > 0) {
   transpose_unsym(A->nrow, A->ncol, A->p, A->i, A->x,
                   AT->nrow, AT->ncol, AT->p, AT->i, AT->x);
   AT->nzmax = A->nzmax;
   ss = new SolverSettings(H, sKKTrhs_init, B, BT, A, AT);
  } else {
   ss = new SolverSettings(H, sKKTrhs_init, B, BT);
  }
#ifdef OPENMP
  ss->ldl_variant = 4;
#else
  ss->ldl_variant = 2;
#endif
  ss->ldl_update_variant = 2;
  ss->solver_mode = 1;
  ss->req_ref_iter = max_iter;
  ss->max_inner_iter = max_iter;
/*
   if (outer_iter_ref > 0 && inner_iter_ref == 0) {
    ss->max_inner_iter = 1; // will be wrong if it is zero
   } else {
    ss->max_inner_iter = inner_iter_ref;
   }
*/

  size_t total_nnz = H->nzmax + B->nzmax;
  double min_v=0, max_v=0;
  //max_min_spmat(H->ncol, H->p, H->i, H->x,max_v,min_v);
  //setting perturbation
  if (auto_reg_en == 2) {
   double H_norm2 = norm_sparse_int(H->ncol, H->p, H->i, H->x, -1, 2);
   auto_perturbation2(H->ncol, total_nnz, H_norm2, max_v - min_v,
                      B->nrow, B->nrow / double(H->ncol), eps_abs,
                      diag_perturb,
                      outer_iter_ref, inner_iter_ref, stop_tol);
  } else if (auto_reg_en == 1) {
   determine_iterations(total_nnz, eps_abs, outer_iter_ref, inner_iter_ref);
   stop_tol = 1e-15;
   double H_norm2 = norm_sparse_int(H->ncol, H->p, H->i, H->x, -1, 2);
   if(H_norm2 < 1)
    diag_perturb = pow(10, -7);
   else
    diag_perturb = pow(10, -9);
  } else{
   outer_iter_ref = inner_iter_ref = max_iter;
  }
  ss->tol_abs = ss->tol_rel = stop_tol;
  ss->req_ref_iter = outer_iter_ref;
  ss->max_inner_iter = inner_iter_ref;
  ss->reg_diag = diag_perturb;
  zero_thresh = diag_perturb;
  //status = analyze_kkt();
  ss->symbolic_analysis();
  delete[]sKKTrhs_init; //TODO: define a new ss construtor
  skkt_col = ss->SM->ncol;
  sKKTrhs = ss->rhs;
  kkt_solution = ss->x;
  descent = kkt_solution;
  lagrange_mult_eq = kkt_solution + H->ncol;
  lagrange_mult = kkt_solution + H->ncol + A->nrow;
  sKKT = ss->A_ord;
#ifdef  REC
  /// compute rec length of constraints.
   //print_csc("B\n",B->ncol,B->p,B->i,B->x);
   compute_recieporical_length(B,b,rec_length);
   //print_vec("rec_length\n",0,B->nrow,rec_length);
#endif
  return status;
 }

 void Nasoq::reset_symbolic_info() { //
 }

 void Nasoq::auto_perturbation(size_t total_nnz, double H_norm2, size_t dim, double &reg_diag, int &outer_iter_ref,
                               int &inner_iter_ref) {
  if (total_nnz < 3e4 || dim < 4100) {
   reg_diag = pow(10, -11);
   inner_iter_ref = 2;
   outer_iter_ref = 2;
   if (H_norm2 < 1) {
    reg_diag = pow(10, -7);
   }
  } else {
   if (H_norm2 < 1) {
    reg_diag = pow(10, -7);
    inner_iter_ref = 1;
    outer_iter_ref = 1;
   } else if (H_norm2 < 1e3) {
    reg_diag = pow(10, -11);
    inner_iter_ref = 9;
    outer_iter_ref = 9;
   } else if (H_norm2 < 5e5) {
    reg_diag = pow(10, -8);
    inner_iter_ref = 0;
    outer_iter_ref = 0;
   } else {
    reg_diag = pow(10, -11);
    inner_iter_ref = 1;
    outer_iter_ref = 1;
   }
  }
 }

 void
 Nasoq::auto_perturbation2(size_t dim, size_t total_nnz, double H_norm2, double range, int n_ineq, double ineq_ratio,
                           double acc_thresh, double &reg_diag, int &outer_iter_ref, int &inner_iter_ref, double &tol) {
  reg_diag = pow(10, -8);
  int iter_low = 2;
  int iter_high = 9;
  int iter_large = 0;
  if (acc_thresh >= 1e-3) {
   iter_low /= 2;
   iter_high /= 2;
  } else if (acc_thresh < 1e-8) {
   iter_low *= 2;
   iter_large *= 2;
   iter_large += 1;
  }
  //std::cout<<acc_thresh<<";"<<iter_low<<";"<<iter_high<<";"<<iter_large<<"\n";
  inner_iter_ref = iter_large;
  outer_iter_ref = iter_large;

  if (total_nnz < 3e5) {
   //std::cout<<"small;";
   inner_iter_ref = iter_low;
   outer_iter_ref = iter_low;
   tol = 1e-10;
//   if (range > 1e3) {
//    //std::cout<<"range;"<<range<<";";
//    diag_perturb = pow(10, -9);
//    tol = 1e-16;
//   } else
   if (H_norm2 < 1) {
    reg_diag = pow(10, -7);
    inner_iter_ref = iter_low;
    outer_iter_ref = iter_low;
    tol = 1e-15;
   } else if (ineq_ratio > 2.5) {
    reg_diag = pow(10, -10);
    tol = 1e-16;
    inner_iter_ref = iter_low;
    outer_iter_ref = iter_low;
   } else if (dim > 1e4 && n_ineq >= 1e4) {
    //std::cout<<"ss;"<<ineq_ratio<<";";
    reg_diag = pow(10, -11);
    inner_iter_ref = iter_high;
    outer_iter_ref = iter_high;
    tol = 1e-15;
   }
  } else {
   if (H_norm2 < 1) {
    reg_diag = pow(10, -7);
    inner_iter_ref = 2;
    outer_iter_ref = 2;
    tol = 1e-15;
   }
  }
 }

 void Nasoq::determine_iterations(size_t total_nnz, double acc_thresh, int &outer_iter_ref, int &inner_iter_ref) {

  if (acc_thresh >= 1e-3) {
   inner_iter_ref = outer_iter_ref = 0;
  } else if (acc_thresh <= 1e-4 && acc_thresh > 1e-7) {
   inner_iter_ref = outer_iter_ref = 1;
  } else if (acc_thresh <= 1e-7 && acc_thresh > 1e-10) {
   inner_iter_ref = outer_iter_ref = 2;
  } else {
   inner_iter_ref = outer_iter_ref = 4;
  }

/*  if(ineq_ratio>2.5 || n_ineq >= 1e4){
  //inner_iter_ref = outer_iter_ref = 3;
  diag_perturb=pow(10,-11);
 }*/

  if (total_nnz <= 20000) {
   inner_iter_ref = outer_iter_ref = 3;
  }
  //std::cout<<acc_thresh<<";"<<iter_low<<";"<<iter_high<<";"<<iter_large<<"\n";

/* if(H_norm2<1){
  diag_perturb=pow(10,-7);
 }*//*else if(ineq_ratio>2.5 || n_ineq >= 1e4){
   diag_perturb = pow(10, -11);
   inner_iter_ref = outer_iter_ref = 9;
  }*/

 }

 void Nasoq::allocate_workspace() {
  // for update func
  //int workspace: 4*super_max + 2*n + 3*supNo
  // double workspace: 2 * super_max*col_max
  size_t upd_int = 4 * ss->max_sup_wid + 2 * H->ncol + 3 * ss->L->nsuper;
  size_t upd_dbl = 2 * ss->max_sup_wid * ss->max_col;

  // workspace for objectives and norms
  size_t nrm_dbl = B->nrow + 4 * H->ncol;


 }

 int Nasoq::set_warm_start(double *x0_in, double *dual_in, const std::vector<int> &actve_set, double *dual_in_eq) {
  //x0 = x0_in;
  for (int i = 0; i < H->ncol; ++i) {
   primal_vars[i] = x0_in[i];
  }
  //setting dual vars
  for (int l = 0; l < A->nrow; ++l) {
   dual_vars_eq[l] = dual_in_eq[l];
  }
  for (int j = 0; j < B->nrow; ++j) {
   dual_vars[j] = dual_in[j];
  }
  n_active = 0;
  for (int k = 0; k < actve_set.size(); ++k) {
   int c_no = actve_set[k];
   active_set.push_back(c_no);
   assert(c_no < B->nrow);
   edit_kkt(1, c_no);
   n_active++;
   used_const[c_no] = 1;
  }
  warm_start = 1;
  return warm_start;
 }

 int Nasoq::set_warm_start_primal(double *x0_in) {
  //x0 = x0_in;
  for (int i = 0; i < H->ncol; ++i) {
   primal_vars[i] = x0_in[i];
  }
  n_active = 0;
  warm_start = 1;
  return warm_start;
 }

 int Nasoq::initialize_x() {
  for (int i = 0; i < H->ncol; ++i) {
   sKKTrhs[i] = -q[i];
  }
  for (int i = H->ncol, j = 0; i < H->ncol + A->nrow; ++i, ++j) {
   sKKTrhs[i] = a[j];
  }
  for (int i = H->ncol + A->nrow; i < sKKT->ncol; ++i) {
   sKKTrhs[i] = .0;
  }
  ss->ldl_update_variant = 2;
  //print_vec("\n RHS: \n",0,sKKT->ncol,sKKTrhs);
  solve_kkt(UPDATE); //unconstrained solution

  ss->ldl_update_variant = 2;
  for (int i = 0; i < H->ncol; ++i) {
   primal_vars[i] = kkt_solution[i];
  }
  for (int i = 0, j = H->ncol; i < A->nrow; ++i, ++j) {
   dual_vars_eq[i] = kkt_solution[j];
  }
/*  for (int k = 15293; k < 16465; ++k) {
  primal_vars[k] = 0;
 }*/
  //print_vec("\nX0: \n",0,H->ncol,primal_vars);
  //print_vec("\nD1: \n",H->ncol,H->ncol+A->nrow,kkt_solution);
  //print_vec("\nD2: \n",H->ncol+A->nrow,H->ncol+A->nrow+B->nrow,kkt_solution);
  return 1;
 }

 int Nasoq::numeric_QP() {
  int nxt_active = -1;
  int nxt_drop = -1;
  int to_add = 1;
  int num_violated = 0;
  double prime_step, dual_step, step;
  solve_type typ_nxt_solve = SOLVE;
  if (warm_start)
   typ_nxt_solve = REFACTOR;
  while (true) {
   if (num_iter > max_iter_nas) {
    is_converged = 0;
    return is_converged;
   }
   if (to_add == 1) {
    num_violated = primal_feasibility(nxt_active);
    //std::cout<<" -> "<< num_violated <<"\n";
    if (num_violated == 0) {
     is_converged = 1;
     return is_converged;
    } else if (nxt_active == -1) {//unbounded!
     std::cout << "All constraints are used\n";
     is_converged = 0;
     return is_converged;
    }
    update_rhs(nxt_active);
    //update_active_set(to_add,nxt_active);
#ifdef VERBOSE
    std::cout<<" Trying with : "<<nxt_active<<" Iter: "<<num_iter<<" ,violated:"<<num_violated<<"\n";
#endif
   }

   /*if(num_iter == 68)
    std::cout<<"-> "<<num_iter<<"\n";*/
/*   int num_neg_d = dual_feasibility();
  assert(num_neg_d == 0);*/
   //std::cout<<"OBJ: "<<compute_primal_obj()<<"\n";
#ifdef CHOLROWMOD
   solve_kkt_cholmod(typ_nxt_solve);
#else
   solve_kkt(typ_nxt_solve);
#endif
   //std::cout<<num_iter<<" ";
//   if(num_iter == 114){
//    int statu;
//    CSC *TMP = ptranspose(ss->A_ord,2,ss->L->IPerm,NULL,0,statu);
//    CSC *TMP2 = ptranspose(TMP,2,NULL,NULL,0,statu);
//    print_csc("\noriginal \n",TMP2->ncol,TMP2->p,TMP2->i,TMP2->x);
//    print_vec("",0,sKKT->ncol,sKKTrhs);
//    print_vec<double >("\nkkt sol\n: ",0,sKKT->ncol,kkt_solution);
//   }
#ifdef VERBOSE
   print_vec<double >("\nlags \n: ",0,B->nrow,lagrange_mult);
    print_vec<double >("\ndual \n: ",0,B->nrow,dual_vars);
    print_vec<double >("\nprimal \n: ",0,H->nrow,primal_vars);
    print_vec<double >("\nkkt sol\n: ",0,sKKT->ncol,kkt_solution);
#endif
/*   if(num_iter == 905)
   printf("\n");
  std::cout<<"--->"<<dual_vars[896]<<"\n";*/
   prime_step = primal_step_length(nxt_active);
   dual_step = dual_step_length(nxt_drop);
   step = prime_step > dual_step ? dual_step : prime_step;
#ifdef VERBOSE
   std::cout<<"Iter: "<<num_iter<<" P-step: "<<prime_step<<" D-step: "
             <<dual_step<<" step:"<<step<<"\n";
#endif
   //std::cout<<"Step sizes: "<<step<<"\n";
   if (step == std::numeric_limits<double>::max()) {
    is_converged = 0;
    return is_converged;
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
#ifdef VERBOSE
    std::cout<<"Iter: "<<num_iter<<" Removed : "<<nxt_drop<<" ,violated:"<<num_violated<<"\n";
#endif
   } else { // Full step
    update_primal(step);// solution = step*descent + solution
    update_dual(step);// dual_vars = step*descent + dual_vars
    dual_vars[nxt_active] = dual_vars[nxt_active] + step;//==-step*-1
    if (step == prime_step) { //update
#ifdef VERBOSE
     std::cout<<"Iter: "<<num_iter<<" Activated : "<<nxt_active<<" ord: "<<
               ss->L->IPerm[nxt_active] <<" ,violated:"<<num_violated<<"\n";
#endif
     to_add = 1;
#ifdef CHOLROWMOD
     edit_kkt_cholmod(to_add, nxt_active);
#else
     edit_kkt(to_add, nxt_active);
#endif

     typ_nxt_solve = UPDATE;
     update_active_set(to_add, nxt_active);
    } else if (step == dual_step) { //downdate
#ifdef VERBOSE
     std::cout<<"Iter: "<<num_iter<<" Removed : "<<nxt_drop<<" ,violated:"<<num_violated<<"\n";
#endif
     to_add = 0;
#ifdef CHOLROWMOD
     edit_kkt_cholmod(to_add, nxt_drop);
#else
     edit_kkt(to_add, nxt_drop);
#endif
     typ_nxt_solve = UPDATE;
     int is_up = update_active_set(to_add, nxt_drop);
     assert(is_up == 1);
#ifdef VERBOSE
     std::cout<<"dropped : "<<nxt_drop<<"\n";
#endif
    }
    assert(n_active == active_set.size());
   }
#ifdef VERBOSE
   print_vec<double >("\nprimal var\n: ",0,H->ncol,primal_vars);
    print_vec<double >("\ndual var\n: ",0,B->nrow,dual_vars);
    std::cout<<"----> end of Iter: "<<num_iter<<"\n";
#endif
   num_iter++;
  }
 }

 void Nasoq::update_dual(double step) {
  for (int i = 0; i < A->nrow; ++i) {
   dual_vars_eq[i] -= step * lagrange_mult_eq[i];
  }
  for (int i = 0; i < active_set.size(); ++i) {
   int c = active_set[i];
   dual_vars[c] -= step * lagrange_mult[c];
  }
  //This is its equivalent
  /*for (int i = 0; i < B->nrow; ++i) {
   dual_vars[i] = dual_vars[i] - step*lagrange_mult[i];
  }
  */
 }

 void Nasoq::update_primal(double step) {
  //TODO using axpy:
  for (int i = 0; i < hessian_size; ++i) {
   primal_vars[i] = primal_vars[i] - step * descent[i];
  }
 }

 int Nasoq::update_active_set(int add_drop, int c_no) {
  std::vector<int>::iterator it;
  if (add_drop) { // add constraint
   n_active++;
   active_set.push_back(c_no);
   used_const[c_no] = 1;
  } else {
   n_active--;
   it = find(active_set.begin(), active_set.end(), c_no);
   if (active_set.size() == 0)
    return 0;
   if (active_set.size() == 1 && active_set[0] == c_no ||
       it != active_set.end() && active_set.size() > 1) {
    active_set.erase(it);
    used_const[c_no] = 0;
    return 1;
   }
   return 0;//c_no not in active_set!!
  }
  return 1;
 }

 int Nasoq::solve() {
  ret_val = nasoq_status::NotConverged;
  int status = 0;
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
   for (int i = 0; i < 4; i++) {
    initialize_x();
    status = numeric_QP();
   }

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
   status = symbolic_QP();
   initialize_x();
#endif
   status = numeric_QP();
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
  compute_objective();

  ret_val = check_solve_status(status);
  detect_solver_name();
  return ret_val;
 }

 nasoq_status Nasoq::check_solve_status(int st) {
  if(st == 0)
   return nasoq_status::Infeasible;
  constraint_sat_norm();
  lagrangian_residual_norm();
  complementarity_norm();
  non_negativity_norm();
  if (lag_res <= eps_abs && cons_sat_norm <= eps_abs &&
      non_negativity_infn <= eps_abs) {//&& complementarity_infn <= eps_abs
   ret_val = nasoq_status::Optimal; //converged
  } else if (cons_sat_norm <= eps_abs) {//low accuracy
   ret_val = nasoq_status::Inaccurate;
  } else{
   ret_val = nasoq_status::NotConverged; //not converged or maybe very inaccurate
  }
  return ret_val;
 }

 void Nasoq::detect_solver_name() {
  if (variant == nasoq_mode::Fixed)
   sol_name = "NASOQ-Fixed";
  else if (variant == nasoq_mode::AUTO)
   sol_name = "NASOQ-AUTO";
  else if (variant == nasoq_mode::Tuned)
   sol_name = "NASOQ-TUNED";
  else
   sol_name = "NASOQ";
 }

 int Nasoq::solve(double *primal_in, double *dual_in, const std::vector<int> &as) {
  int ret_val = 0;
  warm_start = 1;
  qi->tot_st = qi->tic();
  if (scaling > 0) {
   scale_data();
   print_csc("H scaled\n", H->ncol, H->p, H->i, H->x);
  }
  ret_val = symbolic_QP();
  set_warm_start(primal_in, dual_in, as);
  ret_val = numeric_QP();
  if (scaling > 0) {
   unscale_data();
   unscale_solution();
  }
  qi->tot_end = qi->toc();
  qi->tot = qi->elapsed_time(qi->tot_st, qi->tot_end);
  compute_objective();
  constraint_sat_norm();
  lagrangian_residual_norm();
  return ret_val;
 }

 int Nasoq::solve(double *primal_in) {
  int ret_val = 0;
  warm_start = 1;
  qi->tot_st = qi->tic();
  if (scaling > 0) {
   scale_data();
   print_csc("H scaled\n", H->ncol, H->p, H->i, H->x);
  }
  ret_val = symbolic_QP();
  set_warm_start_primal(primal_in);
  ret_val = numeric_QP();
  if (scaling > 0) {
   unscale_data();
   unscale_solution();
  }
  qi->tot_end = qi->toc();
  qi->tot = qi->elapsed_time(qi->tot_st, qi->tot_end);
  compute_objective();
  constraint_sat_norm();
  lagrangian_residual_norm();
  return ret_val;
 }

 int Nasoq::primal_feasibility(int &nxt_add) {
  nxt_add = -1;
  double max_ax = 0;
  int violated_no = 0;
  double *tmp = workspace;
/*  double tmp1,tmp2=10,tmp3=0.1;
 tmp1 = (tmp2= tmp2,tmp3);
 std::cout<<tmp1;*/
  if (B->nrow == 0)
   return violated_no;
  //tmp = B*primal
  //double nnn = norm_sparse_int(B->ncol,B->p,B->i,B->x,0,1);
  //std::cout<<"normx: "<<nnn<<"\n";
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
/*  if(scaling){
  for (int i = 0; i < H->ncol; ++i) {
   tmp[i] *= Einv[i];
  }
 }*/
  for (int i = 0; i < B->nrow; ++i) {
   if (max_ax < std::abs(tmp[i]))
    max_ax = std::abs(tmp[i]);
  }
#ifdef VERBOSE
  print_vec(" Ax: ",0,B->nrow,tmp);
#endif
  //Find the most violent constraint
  double tol = eps_abs;//+ eps_rel*max_ax;
  double max_violated = tol;//||B*x||inf * eabs
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] -= b[i];
   if (tmp[i] > tol) {
    violated_no++;
   }
#ifndef  REC
   if (tmp[i] > max_violated) {
    nxt_add = i;
    max_violated = tmp[i];
   }
#else
   if(tmp[i]*rec_length[i] > max_violated && tmp[i] > tol){
     nxt_add = i;
     max_violated = rec_length[i]*tmp[i];
    }
#endif
  }

#ifdef VERBOSE
  print_vec(" Ax-b: ",0,B->nrow,tmp);
   std::cout<<"max vio "<<tmp[nxt_add]<<"\n";
   std::cout<<"max vio "<<tmp[9]<<"\n";
#endif
  return violated_no;
 }

 int Nasoq::dual_feasibility() {
  int neg_duals = 0;
  for (int i = 0; i < num_active; ++i) {
   if (dual_vars[i] < 0) {
    neg_duals++;
   }
  }
  return neg_duals;
 }

 double Nasoq::dual_step_length(int &nxt_drop) {
  double step = std::numeric_limits<double>::max();
  for (int i = 0; i < B->nrow; ++i) {//TODO: we don't need all range
   if (lagrange_mult[i] > 0 && used_const[i]) {
    double tmp = dual_vars[i] / lagrange_mult[i];
    if (tmp < step) {
     step = tmp;
     nxt_drop = i;
    }
   }
  }
  return step;
 }

 double Nasoq::primal_step_length(int nxt_active) {
  double step = std::numeric_limits<double>::max(), bdotx = 0, ddotb = 0;
  for (int i = BT->p[nxt_active]; i < BT->p[nxt_active + 1]; ++i) {
   bdotx += (BT->x[i] * primal_vars[BT->i[i]]);
   ddotb += (BT->x[i] * descent[BT->i[i]]);
  }
  //std::cout<<"++ : "<<ddotb<<"\n";
  if (ddotb < zero_thresh) {
   return step;
  }
  step = -(b[nxt_active] - bdotx) / ddotb;
  return step;
 }

 int Nasoq::update_rhs(int nxt_active) {
  if (nxt_active < 0 || nxt_active >= B->nrow)
   return 0;
  std::fill_n(sKKTrhs, skkt_col, 0);
  for (int i = BT->p[nxt_active]; i < BT->p[nxt_active + 1]; ++i) {
   sKKTrhs[BT->i[i]] = BT->x[i];
  }
  return 1;
 }

 int Nasoq::solve_kkt_from_scratch() {
  //reseting lagrange multipliers
  for (int i = 0; i < A->nrow; ++i) {
   lagrange_mult_eq[i] = 0;
  }
  for (int i = 0; i < B->nrow; ++i) {
   lagrange_mult[0] = 0;
  }
  double *packed_sol = new double[H->ncol + A->nrow + active_set.size()]();
  //print_vec<double >("\n\nrrhhssss\n: ",0,sKKT->ncol,sKKTrhs);
  build_super_solve_with_eq(H, B, A, sKKTrhs, diag_perturb, active_set, packed_sol,
                            outer_iter_ref, inner_iter_ref, stop_tol);

/*    build_super_solve_with_eq_mkl(H,B,A,sKKTrhs,diag_perturb,active_set,packed_sol,
   outer_iter_ref,inner_iter_ref,stop_tol);*/
  for (int i = 0; i < H->ncol; ++i) {
   descent[i] = packed_sol[i];
  }
  for (int i = 0; i < A->nrow; ++i) {
   lagrange_mult_eq[i] = packed_sol[H->ncol + i];
  }
  for (int i = 0; i < active_set.size(); ++i) {
   lagrange_mult[active_set[i]] = packed_sol[H->ncol + A->nrow + i];
  }
  delete[]packed_sol;
  return 1;
 }

 int Nasoq::solve_kkt(solve_type s_type) {
  int ret_val = 0;
#if 0 //ndef SYM_REMOV
  solve_kkt_from_scratch();
#else
  qi->fct_st = qi->tic();
  if (s_type == SOLVE) {
   ss->solve_only();
  } else if (s_type == UPDATE) {
/*   if(qi->sw){ // For measuring the effect of initial factorization
   qi->ini_st = qi->tic();
  }*/
   ss->update_factorization();
/*   if(qi->sw){
   qi->ini_end = qi->toc();
   qi->init = qi->elapsed_time(qi->ini_st,qi->ini_end);
   qi->sw = false;
  }*/
   //ss->check_ldlt_factor();
   ss->solve_only(); //TODO replace with update_solve
  } else { // refactor
   ss->numerical_factorization();
   //ss->check_ldlt_factor();
   ss->solve_only();
  }
  qi->fct_end = qi->toc();
  qi->factt += qi->elapsed_time(qi->fct_st, qi->fct_end);
#endif
/*  double bwd_err = ss->backward_error();
 std::cout<<bwd_err<<"\n";*/
  return ret_val;
 }

 void Nasoq::edit_kkt(int add_del, int nxt_col) {
  std::vector<int> sn_list;
  sn_list.push_back(nxt_col);


  ss->add_del_matrix_qp(add_del, sn_list);
 }

 double Nasoq::compute_primal_obj() {
  primal_obj = 0;
  int status = 0;
  double alp[2] = {1.0, 0};
  double bet[2] = {0.0, 0};
  //double *tmp = workspace;
  double *tmp = new double[H->ncol]();
  spmv_csc_sym_one_int(H->ncol, H->p, H->i, H->x, -1, alp, bet,
                       1, primal_vars, tmp);
  /*CSC *HTT = ptranspose(H,2,NULL,NULL,0,status);
  primal_obj = quad_form(HTT,primal_vars);*/
  primal_obj = 0.5 * dot(H->ncol, tmp, primal_vars);
  primal_obj += dot(H->ncol, q, primal_vars);
  delete[]tmp;
  return primal_obj;
 }

 double Nasoq::quad_form(const CSC *P, const double *x) {
  double quad_form = 0.;
  int i, j, ptr;                                // Pointers to iterate over
  // matrix: (i,j) a element
  // pointer

  for (j = 0; j < P->ncol; j++) {                      // Iterate over columns
   for (ptr = P->p[j]; ptr < P->p[j + 1]; ptr++) { // Iterate over rows
    i = P->i[ptr];                                // Row index

    if (i == j) {                                 // Diagonal element
     quad_form += (double) .5 * P->x[ptr] * x[i] * x[i];
    } else if (i < j) {                             // Off-diagonal element
     quad_form += P->x[ptr] * x[i] * x[j];
    } else {                                        // Element in lower diagonal
     // part
#ifdef PRINTING
     c_eprint("quad_form matrix is not upper triangular");
#endif /* ifdef PRINTING */
     return -1;
    }
   }
  }
  return quad_form;
 }

 double Nasoq::compute_dual_obj() {
  dual_obj = 0;
  double *tmp = workspace;
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] -= b[i];
  }
  dual_obj = dot(B->nrow, tmp, dual_vars);
  return dual_obj;
 }

 double Nasoq::complementarity_norm() {
  complementarity_infn = 0;
  double *tmp = workspace;
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] = dual_vars[i] * (tmp[i] - b[i]);
  }
  complementarity_infn = norm_dense(B->nrow, 1, tmp, 0);

/*  if(A->nrow>0){
  spmv_csc_small(A->nrow,A->ncol,A->p, A->i, A->x, primal_vars, tmp);
  for (int i = 0; i < A->nrow; ++i) {
   tmp[i]= dual_vars_eq[i] * (tmp[i]-a[i]);
  }
  double constraint_sat_norm_eq = norm_dense(A->nrow,1,tmp,0);
  complementarity_infn = std::max(complementarity_infn, constraint_sat_norm_eq);
 }*/
  return complementarity_infn;
 }

 double Nasoq::compute_objective() {
  objective = compute_primal_obj() + compute_dual_obj();
  return objective;
 }

 double Nasoq::lagrangian_residual_norm() {
  lag_res = 0;
  double alp[2] = {1.0, 0};
  double bet[2] = {0.0, 0};
  double *Hx = workspace;
  double *B_eq_dual = workspace + H->ncol;
  double *Bdual = workspace + 2 * H->ncol;
  double *res_norm = workspace + 3 * H->ncol;
  spmv_csc_sym_one_int(H->ncol, H->p, H->i, H->x, -1, alp, bet,
                       1, primal_vars, Hx);
  //print_vec("lag ",0,B->ncol,Hx);
  if (A->nrow > 0) {
   spmv_csc_small(AT->nrow, AT->ncol, AT->p, AT->i, AT->x,
                  dual_vars_eq, B_eq_dual);
  } else {
   for (int i = 0; i < H->ncol; ++i) {
    B_eq_dual[i] = 0;
   }
  }
  //print_vec("Beq_dual",0,H->ncol,B_eq_dual);
  //double nneq = norm_dense(H->ncol,1,B_eq_dual,0);
  if (B->nrow > 0) {
   spmv_csc_small(BT->nrow, BT->ncol, BT->p, BT->i, BT->x,
                  dual_vars, Bdual);
  } else {
   for (int i = 0; i < H->ncol; ++i) {
    Bdual[i] = 0;
   }
  }
  //double nnieq = norm_dense(H->ncol,1,Bdual,0);
  //print_vec("Bdual",0,H->ncol,Bdual);
  for (int i = 0; i < H->ncol; ++i) {
   res_norm[i] = Hx[i] + q[i] + Bdual[i] + B_eq_dual[i];
  }
  //print_vec("res_norm",0,H->ncol,res_norm);
  lag_res = norm_dense(H->ncol, 1, res_norm, 0);
  //std::cout<<";"<<nneq<<";"<<nnieq<<";"<<lag_res<<"\n";
  return lag_res;
 }

 double Nasoq::constraint_sat_norm() {
  cons_sat_norm = 0;
  double *tmp = workspace;
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
  for (int i = 0; i < B->nrow; ++i) {
   double diff = tmp[i] - b[i];
   //tmp[i] = std::min(tmp[i]-b[i],0.0);
   //tmp[i] = std::max(b[i] - tmp[i], 0.0); //Equal to top
   tmp[i] = diff < 0.0 ? 0.0 : diff;
//   if(diff > 0){
//    std::cout<<"dddddd: "<<i<<" ->" <<diff<<"\n";
//   }

  }
  cons_sat_norm = norm_dense(B->nrow, 1, tmp, 0);

  if (A->nrow > 0) {
   spmv_csc_small(A->nrow, A->ncol, A->p, A->i, A->x, primal_vars, tmp);
   for (int i = 0; i < A->nrow; ++i) {
    tmp[i] = tmp[i] - a[i];
    //tmp[i] = std::min(tmp[i]-b[i],0.0);
    //tmp[i] = std::max(b[i] - tmp[i], 0.0); //Equal to top
    //tmp[i] = diff < 0.0 ? 0.0 : diff;
   }
   double constraint_sat_norm_eq = norm_dense(A->nrow, 1, tmp, 0);
   cons_sat_norm = std::max(cons_sat_norm, constraint_sat_norm_eq);
  }
  return cons_sat_norm;
 }

 double Nasoq::dual_FB_norm() {
  dual_FB_a = new double[B->ncol];
  dual_FB_b = new double[B->ncol];
  dual_FB = new double[B->ncol]();
  double *Bdual;
  double *HdivBdual;
  double *Hdivq;
  if (H->ncol >= B->nrow) {
   Bdual = workspace;
   HdivBdual = workspace + H->ncol;
   Hdivq = workspace + 2 * H->ncol;
  } else {
   Bdual = workspace;
   HdivBdual = new double[B->nrow];
   Hdivq = new double[B->nrow];
  }
  spmv_csc_small(BT->nrow, BT->ncol, BT->p, BT->i, BT->x,
                 dual_vars, Bdual); //Bdual = B^T*dual_vars
  SolverSettings *s_fb = new SolverSettings(H, Bdual);
  s_fb->symbolic_analysis();
  s_fb->numerical_factorization();
  s_fb->solve_only(); //x = H^{-1}*Bdual
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x,
                 s_fb->x, Hdivq); //Hdivq = B*x
  //print_vec("B*dual: ",0,B->nrow,Hdivq);
  for (int j = 0; j < H->ncol; ++j) {
   s_fb->rhs[j] = q[j];
  }
  s_fb->solve_only(); // x = H^{-1}*q
  //print_vec("B*dual: ",0,B->ncol,s_fb->x);
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x,
                 s_fb->x, HdivBdual); //HdivBdual = B*x
  //print_vec("B*dual: ",0,B->nrow,HdivBdual);
  for (int i = 0; i < B->nrow; ++i) {
   dual_FB_a[i] = dual_vars[i];
   dual_FB_b[i] = Hdivq[i] + HdivBdual[i] + b[i];
  }
  //print_vec("B*dual: ",0,B->nrow,dual_FB_a);
  //print_vec("B*dual: ",0,B->nrow,dual_FB_b);
  //print_vec("b: ",0,)
  Fischer_Burmeister_func(B->nrow, dual_FB_a, dual_FB_b, dual_FB);
  //print_vec("B*dual: ",0,B->nrow,dual_FB);
  double fb_norm = norm_dense(B->nrow, 1, dual_FB, 0);
  delete s_fb;
  delete[] dual_FB;
  delete[] dual_FB_a;
  delete[] dual_FB_b;
  if (H->ncol < B->nrow) {
   delete[]HdivBdual;
   delete[]Hdivq;
  }
  return fb_norm;
 }

 double Nasoq::primal_FB_norm() {
  double *tmp = workspace;
  primal_FB_a = new double[B->nrow];
  primal_FB_b = new double[B->nrow];
  primal_FB = new double[B->nrow]();
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
  for (int i = 0; i < B->nrow; ++i) {
   primal_FB_a[i] = dual_vars[i];
   primal_FB_b[i] = -tmp[i] + b[i];
  }
  Fischer_Burmeister_func(B->nrow, primal_FB_a, primal_FB_b, primal_FB);
  double fb_norm = norm_dense(B->nrow, 1, primal_FB, 0);
  delete[]primal_FB;
  delete[]primal_FB_a;
  delete[]primal_FB_b;
  return fb_norm;
 }

 double Nasoq::non_negativity_norm() {
  double *tmp = workspace;
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] = dual_vars[i] > 0.0 ? 0 : dual_vars[i];
  }
  non_negativity_infn = norm_dense(B->nrow, 1, tmp, 0);
  return non_negativity_infn;
 }

 void Nasoq::print(int vb_level) {
  if (vb_level >= 1) {
   print_vec<double>("\nprimal var\n: ", 0, H->ncol, primal_vars);
   print_vec<double>("\ndual var\n: ", 0, B->nrow, dual_vars);
   print_vec<double>("\ndual equality var\n: ", 0, A->nrow, dual_vars_eq);
  }
  if (is_converged)
   std::cout << "\nConverged in " << num_iter << " .";
  std::cout << "\n objective is (primal, dual, total):" << primal_obj << ", " << dual_obj << ", " << objective << "\n";
 }

 void Nasoq::export_to_file(std::string p_name, std::string vec_head) {
  std::string primal_name = p_name + "_primal.txt";
  write_vector(primal_name, H->ncol, primal_vars, vec_head);
  std::string dual_name = p_name + "_dual.txt";
  write_vector(dual_name, B->nrow, dual_vars);
  double info[5] = {};
  info[0] = num_iter;
  info[1] = primal_obj;
  info[2] = dual_obj;
  info[3] = objective;
  info[4] = qi->tot;
  std::string other_name = p_name + "_info.txt";
  write_vector(other_name, 5, info);
 }

 void Nasoq::print_log() {

  std::cout << eps_abs << "," << outer_iter_ref << "," << inner_iter_ref << ",";
  std::cout << stop_tol << "," << diag_perturb << "," << ret_val << "," << num_iter << ",";
  std::cout << qi->tot << "," << active_set.size() << "," << cons_sat_norm << ",";
  std::cout << lag_res << "," << primal_obj << "," << dual_obj << "," << objective << ",";
  std::cout << non_negativity_infn << "," << complementarity_infn << ",";
  //std::cout<<qi->init<<",";
 }

 void Nasoq::print_qp_range(int verbose) {
  double max_v = 0, min_v = 0, avg = 0, var = 0;
  double max_max = -std::numeric_limits<double>::max(),
    min_min = std::numeric_limits<double>::max(),
    range = 0;
  max_min_vector(H->nzmax, H->x, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  min_min = std::min(min_min, min_v);
  max_max = std::max(max_max, max_v);
  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;

  max_min_vector(H->ncol, q, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  min_min = std::min(min_min, min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(H->ncol, primal_vars, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  //min_min = std::min(min_min,min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(A->nzmax, A->x, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  min_min = std::min(min_min, min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(A->nrow, a, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  min_min = std::min(min_min, min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(A->nrow, dual_vars_eq, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  //min_min = std::min(min_min,min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(B->nzmax, B->x, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  min_min = std::min(min_min, min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(B->nrow, b, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  min_min = std::min(min_min, min_v);
  max_max = std::max(max_max, max_v);

  max_v = 0;
  min_v = 0;
  avg = 0;
  var = 0;
  max_min_vector(B->nrow, dual_vars, max_v, min_v, avg, var);
  if (verbose)
   std::cout << max_v << "," << min_v << "," << avg << "," << var << ",";
  //min_min = std::min(min_min,min_v);
  max_max = std::max(max_max, max_v);
  range = max_max / min_min;
  std::cout << min_min << "," << max_max << "," << range << ",";
 }

 int Nasoq::scale_data_2() {
  const double max_scale = 1e4;
  double min_h, max_h, avg_h, var_h;
  max_min_sparse_matrix(H, max_h, min_h, avg_h, var_h);
  inv_qp_scalar = (std::min(sqrt(max_h), max_scale));
  inv_qp_scalar = 1;
  qp_scalar = 1.0 / inv_qp_scalar;
  scale_sparse_matrix(H, qp_scalar * qp_scalar);
  scale_vector(H->ncol, q, qp_scalar);
  scale_sparse_matrix(A, qp_scalar);
  scale_sparse_matrix(B, qp_scalar);

  return 1;
 }

 int Nasoq::unscale_data_2() {
  scale_sparse_matrix(H, inv_qp_scalar * inv_qp_scalar);
  scale_vector(H->ncol, q, inv_qp_scalar);
  scale_sparse_matrix(A, inv_qp_scalar);
  scale_sparse_matrix(AT, inv_qp_scalar);
  scale_sparse_matrix(B, inv_qp_scalar);
  scale_sparse_matrix(BT, inv_qp_scalar);

  return 1;
 }

 int Nasoq::unscale_solution_2() {
  scale_vector(H->ncol, primal_vars, qp_scalar);
  //scale_vector(A->nrow,dual_vars_eq,qp_scalar);
  //scale_vector(B->nrow,dual_vars,qp_scalar);
  return 1;
 }

 int Nasoq::scale_data() {
  // Scale KKT matrix
  //
  //    [ P   A'  B']
  //    [ A   0   0 ]
  //    [ B   0   0 ]
  // with diagonal matrix
  //
  //  S = [ D       ]
  //      [    E    ]
  //      [       F ]

  int i;          // Iterations index
  int n, me, m;       // Number of constraints and variables
  double c_temp;     // Cost function scaling
  double inf_norm_q; // Infinity norm of q

  n = H->ncol;
  me = A->nrow;
  m = B->nrow;
  D = new double[n];
  E = new double[me];
  F = new double[m];
  Dinv = new double[n];
  Einv = new double[me];
  Finv = new double[m];
  double *D_temp = workspace; //2n + m + me
  double *D_temp_A = workspace + n;
  double *E_temp = workspace + 2 * n;
  double *F_temp = workspace + 2 * n + me;

  // Initialize scaling to 1
  c = 1.0;
  vec_set_scalar(D, 1., n);
  vec_set_scalar(Dinv, 1., n);
  vec_set_scalar(E, 1., me);
  vec_set_scalar(Einv, 1., me);
  vec_set_scalar(F, 1., m);
  vec_set_scalar(Finv, 1., m);

  for (i = 0; i < scaling; i++) {
   //
   // First Ruiz step
   //

   // Compute norm of KKT columns
   compute_inf_norm_cols_KKT(H, A, B,
                             D_temp, D_temp_A,
                             E_temp, F_temp, n);
   //print_vec("vect:\n",0,n,D_temp);
   // Set to 1 values with 0 norms (avoid crazy scaling)
   limit_scaling(D_temp, n);
   limit_scaling(E_temp, me);
   limit_scaling(F_temp, m);

   // Take square root of norms
   vec_ew_sqrt(D_temp, n);
   vec_ew_sqrt(E_temp, me);
   vec_ew_sqrt(F_temp, m);

   // Divide scalings D and E by themselves
   vec_ew_recipr(D_temp, D_temp, n);
   vec_ew_recipr(E_temp, E_temp, me);
   vec_ew_recipr(F_temp, F_temp, m);

   // Equilibrate matrices P and A and vector q
   // P <- DPD
   mat_premult_diag(H, D_temp);
   mat_postmult_diag(H, D_temp);

   // A <- EAD
   mat_premult_diag(A, E_temp);
   mat_postmult_diag(A, D_temp);

   // B <- FBD
   mat_premult_diag(B, F_temp);
   mat_postmult_diag(B, D_temp);

   // q <- Dq
   vec_ew_prod(D_temp, q, q, n);

   // Update equilibration matrices D and E
   vec_ew_prod(D, D_temp, D, n);
   vec_ew_prod(E, E_temp, E, me);
   vec_ew_prod(F, F_temp, F, m);

   //
   // Cost normalization step
   //

   // Compute avg norm of cols of P
   mat_inf_norm_cols_sym_triu(H, D_temp);
   c_temp = vec_mean(D_temp, n);

   // Compute inf norm of q
   inf_norm_q = vec_norm_inf(q, n);

   // If norm_q == 0, set it to 1 (ignore it in the scaling)
   // NB: Using the same function as with vectors here
   limit_scaling(&inf_norm_q, 1);

   // Compute max between avg norm of cols of P and inf norm of q
   c_temp = std::max(c_temp, inf_norm_q);

   // Limit scaling (use same function as with vectors)
   limit_scaling(&c_temp, 1);

   // Invert scaling c = 1 / cost_measure
   c_temp = 1. / c_temp;

   // Scale P
   mat_mult_scalar(H, c_temp);

   // Scale q
   vec_mult_scalar(q, c_temp, n);

   // Update cost scaling
   c *= c_temp;
  }
  // Store cinv, Dinv, Einv
  cinv = 1. / c;
  vec_ew_recipr(D, Dinv, n);
  vec_ew_recipr(E, Einv, me);
  vec_ew_recipr(F, Finv, m);

  // Scale problem vectors l, u
  vec_ew_prod(E, a, a, me);
  vec_ew_prod(F, b, b, m);
  //vec_ew_prod(E, work->data->u, work->data->u, work->data->m);

  return 0;
 }

 int Nasoq::unscale_data() {
  // Unscale cost
  mat_mult_scalar(H, cinv);
  mat_premult_diag(H, Dinv);
  mat_postmult_diag(H, Dinv);
  vec_mult_scalar(q, cinv, H->ncol);
  vec_ew_prod(Dinv, q, q, H->ncol);

  // Unscale constraints
  if (A->nrow > 0) {
   mat_premult_diag(A, Einv);
   mat_postmult_diag(A, Dinv);

   mat_premult_diag(AT, Dinv);
   mat_postmult_diag(AT, Einv);

   vec_ew_prod(Einv, a, a, A->nrow);
  }

  mat_premult_diag(B, Finv);
  mat_postmult_diag(B, Dinv);

  mat_premult_diag(BT, Dinv);
  mat_postmult_diag(BT, Finv);

  vec_ew_prod(Finv, b, b, B->nrow);

  //vec_ew_prod(Einv, work->data->u, work->data->u, work->data->m);

  return 0;
 }

 int Nasoq::unscale_solution() {
  // primal
  vec_ew_prod(D,
              primal_vars,
              primal_vars,
              H->ncol);

  // dual
  vec_ew_prod(E,
              dual_vars_eq,
              dual_vars_eq,
              A->nrow);
  vec_mult_scalar(dual_vars_eq, cinv, A->nrow);

  vec_ew_prod(F,
              dual_vars,
              dual_vars,
              B->nrow);
  vec_mult_scalar(dual_vars, cinv, B->nrow);

  return 0;
 }

 int Nasoq::numeric_QP_batch() {
  int nxt_active = -1;
  std::vector<idx_val> nxt_list;
  int nxt_drop = -1;
  int to_add = 1;
  int num_violated = 0;
  double prime_step, dual_step, step;
  solve_type typ_nxt_solve = SOLVE;
  while (true) { // Outer iterations
   if (num_iter > max_iter_nas) {
    is_converged = 0;
    return is_converged;
   }

   //Add all violent constraints up to batch size
   if (to_add == 1) {
    num_violated = primal_feasibility_batch(nxt_list);
    if (num_violated == 0) {
     compute_objective();
     is_converged = 1;
     return is_converged;
    } else if (nxt_list.size() == 0) {//unbounded!
     std::cout << "All constraints are used\n";
     is_converged = 0;
     return is_converged;
    }
    edit_kkt_batch(to_add, nxt_list);
    nxt_active = nxt_list[nxt_list.size() - 1].idx;
    update_rhs(nxt_active);
    typ_nxt_solve = UPDATE;
    //solve_kkt(typ_nxt_solve);
    //update_active_set_batch(to_add,nxt_list);
    to_add = 0;
   }
   // Iterartes till no negative multiplier
   while (!to_add) { //Inner iterations
    int neg_dual_idx = -1;
    solve_kkt(typ_nxt_solve);
    int num_neg = dual_feasibility_inner(neg_dual_idx);
    if (num_neg == 0) {
     to_add = 0;
     break;
    }
    edit_kkt(to_add, neg_dual_idx);
    typ_nxt_solve = UPDATE;
    //update the batch list
    remove_from_list(nxt_list, neg_dual_idx);
   }
#ifdef VERBOSE
   print_vec<double >("\nkkt sol\n: ",0,sKKT->ncol,kkt_solution);
#endif
   //here we can add the constraints to active-set
   update_active_set_batch(0, nxt_list);
   prime_step = primal_step_length_batch(nxt_list);
   dual_step = dual_step_length(nxt_drop);
   step = prime_step > dual_step ? dual_step : prime_step;
#ifdef VERBOSE
   std::cout<<"P-step: "<<prime_step<<" D-step: "
             <<dual_step<<" step:"<<step<<"\n";
#endif
   if (step == std::numeric_limits<double>::max()) {
    is_converged = 0;
    return is_converged;
   } else if (prime_step == std::numeric_limits<double>::max()
              && step == dual_step) { // Dual step
    //update_dual(step);
    to_add = 0;
    edit_kkt(to_add, nxt_drop);
    typ_nxt_solve = UPDATE;
    int is_up = update_active_set(to_add, nxt_drop);
    assert(is_up == 1);
   } else { // Full step
    update_primal(step);// solution = step*descent + solution
    update_dual(step);// dual_vars = step*descent + dual_vars
    dual_vars[nxt_active] = dual_vars[nxt_active] + step;//==-step*-1
    if (step == prime_step) { //update
     to_add = 1;
     edit_kkt(to_add, nxt_active);
     //typ_nxt_solve = UPDATE;
     update_active_set(to_add, nxt_active);
    } else if (step == dual_step) { //downdate
     to_add = 0;
     edit_kkt(to_add, nxt_drop);
     typ_nxt_solve = UPDATE;
     int is_up = update_active_set(to_add, nxt_drop);
     assert(is_up == 1);
    }
    assert(n_active == active_set.size());
   }
   //num_iter++;
   //}
#ifdef VERBOSE
   print_vec<double >("\nprimal var\n: ",0,H->ncol,primal_vars);
    print_vec<double >("\ndual var\n: ",0,B->nrow,dual_vars);
    std::cout<<"----> end of Iter: "<<num_iter<<"\n";
#endif
   num_iter++;
  }
 }

 int Nasoq::edit_kkt_batch(int add_del, std::vector<idx_val> nxt_list) {
  std::vector<int> sn_list;
  for (int i = 0; i < nxt_list.size() - 1; ++i) {
   sn_list.push_back(nxt_list[i].idx);
  }
  ss->add_del_matrix_qp(add_del, sn_list);
  return 1;
 }

 double Nasoq::primal_step_length_batch(const std::vector<idx_val> &nxt_list) {
  double step = std::numeric_limits<double>::max(), bdotx = 0, ddotb = 0;
  double cur_step;
  int nxt_active = -1;
  for (int j = 0; j < nxt_list.size(); ++j) {
   nxt_active = nxt_list[j].idx;
   bdotx = 0;
   ddotb = 0;
   for (int i = BT->p[nxt_active]; i < BT->p[nxt_active + 1]; ++i) {
    bdotx += (BT->x[i] * primal_vars[BT->i[i]]);
    ddotb += (BT->x[i] * descent[BT->i[i]]);
   }
   if (ddotb < zero_thresh) {
    return std::numeric_limits<double>::max();
   }
   cur_step = -(b[nxt_active] - bdotx) / ddotb;
   if (cur_step < step) {
    step = cur_step;
   }
  }
  return step;
 }

 int Nasoq::primal_feasibility_batch(std::vector<idx_val> &pair_list) {
  pair_list.clear();
  int violated_no = 0;
  double *tmp = workspace;
  double max_violated = zero_thresh;
  //tmp = B*primal
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
#ifdef VERBOSE
  print_vec(" Ax: ",0,B->nrow,tmp);
#endif
  //Find the violent constraints
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] -= b[i];
   if (tmp[i] > zero_thresh) {
    violated_no++;
    idx_val tmp_pair(i, tmp[i]);
    pair_list.push_back(tmp_pair);
   }
  }
  if (violated_no == 0)
   return violated_no;
  std::sort(pair_list.begin(), pair_list.end(), CMP);
  return violated_no;
 }

 int Nasoq::dual_feasibility_inner(int &most_neg) {
  int neg_duals = 0;
  most_neg = -1;
  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < B->nrow; ++i) {
   if (lagrange_mult[i] < 0) {//FIXME zero_threshold?
    neg_duals++;
    if (lagrange_mult[i] < min) {
     min = lagrange_mult[i];
     most_neg = i;
    }
   }
  }
  return neg_duals;
 }

 int Nasoq::update_active_set_batch(int add_drop, const std::vector<idx_val> nxt_lst) {
  std::vector<int>::iterator it;
  if (add_drop) { // add constraint
   for (int i = 0; i < nxt_lst.size() - 1; ++i) {
    n_active++;
    active_set.push_back(nxt_lst[i].idx);
    used_const[nxt_lst[i].idx] = 1;
   }
  } else {
   int c_no = 0; //FIXME: else part is wrong
/*   n_active--;
  it = find (active_set.begin(), active_set.end(), c_no);
  if(active_set.size()==0)
   return 0;
  if(active_set.size() == 1 && active_set[0] == c_no ||
     it != active_set.end() && active_set.size()>1){
   active_set.erase(it);
   used_const[c_no] = 0;
   return 1;
  }*/
   return 0;//c_no not in active_set!!
  }
  return 1;
 }
}
