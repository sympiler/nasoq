//
// Created by Shujian Qian on 2020-10-30.
//

#include "nasoq/QP/linear_solver_wrapper.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <nasoq/common/Sym_BLAS.h>

#include "nasoq/common/Norm.h"
#include "nasoq/common/Transpose.h"
#include "nasoq/common/TreeUtils.h"
#include "nasoq/common/Util.h"
#include "nasoq/gmres/mgmres.hpp"
#include "nasoq/ldl/ldlt_check.h"
#include "nasoq/ldl/parallel_blocked_ldlt.h"
#include "nasoq/ldl/parallel_blocked_ldlt_02.h"
#include "nasoq/ldl/parallel_blocked_ldlt_03.h"
#include "nasoq/ldl/Parallel_simplicial_ldl.h"
#include "nasoq/ldl/Parallel_update_ldl_02_2.h"
#include "nasoq/ldl/Parallel_update_simplicial.h"
#include "nasoq/ldl/Serial_blocked_ldl.h"
#include "nasoq/ldl/Serial_blocked_ldl_02_2.h"
#include "nasoq/ldl/serial_simplicial_ldl.h"
#include "nasoq/ldl/Serial_update_ldl.h"
#include "nasoq/ldl/Serial_update_ldl_static.h"
#include "nasoq/ldl/Serial_update_simplicial_ldl.h"
#include "nasoq/linear_solver/solve_phase.h"
#include "nasoq/symbolic/symbolic_phase.h"



namespace nasoq {

 profiling_solver_info::profiling_solver_info(int nt) : fact_time(0), analysis_time(0),
                                                        solve_time(0), iter_time(0),
                                                        ordering_time(0), update_time(0),
                                                        piv_reord(0) {
  timing_chol = new double[4 + nt]();
 }

 profiling_solver_info::profiling_solver_info(int nt, double ft, double at, double st, double it, double ot, double ut,
                                              double pr) :
   fact_time(ft), analysis_time(at),
   solve_time(st), iter_time(it),
   ordering_time(ot), update_time(ut),
   piv_reord(pr) {
  timing_chol = new double[4 + nt]();
 }

 profiling_solver_info::~profiling_solver_info() {
  delete[]timing_chol;
 }

 std::chrono::time_point<std::chrono::system_clock> profiling_solver_info::tic() {
  return std::chrono::system_clock::now();
 }

 std::chrono::time_point<std::chrono::system_clock> profiling_solver_info::toc() {
  return std::chrono::system_clock::now();
 }

 double profiling_solver_info::elapsed_time(std::chrono::time_point<std::chrono::system_clock> beg,
                                            std::chrono::time_point<std::chrono::system_clock> lst) {
  double ret = 0;
  elapsed_seconds = lst - beg;
  ret = elapsed_seconds.count();
  return ret;
 }

 void profiling_solver_info::print_profiling() {
  std::cout << "analysis time: " << analysis_time << ";";
  std::cout << "fact time: " << fact_time << ";";
  std::cout << "update time: " << update_time << ";";
  std::cout << "reordering pivot time: " << piv_reord << ";";
  std::cout << "solve time: " << solve_time << ";";
 }

 SolverSettings::SolverSettings(CSC *Amat, double *rhs_in) {
  default_setting();
  A = Amat;
  rhs = rhs_in;
  A_ord = NULL;//new CSC;
  AT_ord = NULL; //new CSC;
  SM = new CSC;
  psi = new profiling_solver_info(num_thread);

  solver_mode = 0;//basic mode
  remove_trans = 0;
  x = NULL;
  B = NULL;
  C = NULL;
  a_consistent = 1;
  to_del = 0;
  visible_sn = NULL;
 }

 SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat) :
   SolverSettings(Amat, rhs_in) { //This an expert routine
  B = Bmat;
  BT = new CSC;
  transpose_unsym(B->nrow, B->ncol, B->p, B->i, B->x,
                  BT->nrow, BT->ncol, BT->p, BT->i, BT->x);
  BT->stype = B->stype;
  BT->xtype = B->xtype;
  BT->packed = B->packed;
  BT->nz = B->nz;
  BT->sorted = B->sorted;
  remove_trans = 1;
  solver_mode = 1;
 }

 SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, double *b) :
   SolverSettings(Amat, rhs_in, Bmat) { //This an expert routine
  extra_rhs = b;
 }

 SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, CSC *BTmat) :
   SolverSettings(Amat, rhs_in) { //This an expert routine
  B = Bmat;
  BT = BTmat;
  solver_mode = 1;
 }

 SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, CSC *BTmat, CSC *Cmat, CSC *CTmat) :
   SolverSettings(Amat, rhs_in, Bmat, BTmat) { //This an expert routine
  C = Cmat;
  CT = CTmat;
  solver_mode = 1;
 }

 SolverSettings::~SolverSettings() {
  allocateAC(A_ord, 0, 0, 0, FALSE);
  allocateAC(AT_ord, 0, 0, 0, FALSE);
  allocateLC(L, FALSE);
  if (simplicial_alloc) {
   delete[]L->p_s;
   delete[]L->x_s;
  }
  delete psi;
  delete[]ws;
  delete[]ws_int;
  delete[]ws_zeroed;
  //delete []pinv;
  delete[]perm_piv;
  delete[]marked;
  delete[]valL;
  delete[]d_val;
  delete[]visible_cnt;

  if (remove_trans) {
   allocateAC(BT, 0, 0, 0, FALSE);
  }

  if (solver_mode == 1 || solver_mode == 2) {
#ifdef SYM_REMOV
   delete []child_sn_ptr;
    delete []child_sn_no;
    delete []num_sn_child;
#endif
   delete[]atree;
   delete[]visible_sn;
#if 0
   delete []child_ptr;
    delete []child_no;
    delete []num_child;
    delete []etree_mod;
#endif
  }

  delete[]level_ptr;
  delete[]par_ptr;
  delete[]par_set;
  if (simplicial_alloc) {
   //delete []level_ptr_s;
   delete[]par_ptr_s;
   delete[]par_set_s;
  }

  if (solver_mode == 2) {
   delete[]sm_rhs;
  }
  if (solver_mode == 1) {
   delete[]sm_rhs;
   delete[]sm_solution;
  }
  delete[]extra_cols;

  if (num_thread > thread_thresh && solver_mode == 1) {
   delete[]s_level_ptr;
   delete[]s_level_set;
  }

  // deleting matrices
  delete[]L->ColCount;
  delete L;
  if (solver_mode != 0)
   allocateAC(SM, 0, 0, 0, FALSE);
  else
   delete SM;
 }

 void SolverSettings::default_setting() {
#ifdef OPENMP
  ldl_variant = 4;
#else
  ldl_variant = 2;
#endif
  ldl_update_variant = 2;
#if defined(MKL_BLAS)
  num_thread = mkl_get_max_threads();
  MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
#elif defined(OPENBLAS)
  num_thread = openblas_get_num_procs();
  openblas_set_num_threads(1);
#else
#error couldn't determine BLAS implementation
#endif
  chunk = 1;
  cost_param = num_thread;
  level_param = -3;
  final_seq_node = 50;// should be greater than 1
  n_relax[0] = 4;
  n_relax[1] = 16;
  n_relax[2] = 48;
  z_relax[0] = 0.8;
  z_relax[1] = 0.1;
  z_relax[2] = 0.05;
  //z_relax[0]=0.9; z_relax[1]=0.5; z_relax[2]=0.05;
  //refinemen
  req_ref_iter = 2;
  max_iter = req_ref_iter;
  max_inner_iter = 2;
  tol_abs = 1e-15;
  tol_rel = 1e-15;
  //
  reg_diag = 1e-8;
  is_super = 1;
  regularization = 1;
  thread_thresh = 25;
  simplicial_alloc = 0;
 }

 int SolverSettings::build_super_matrix() {
  int status=0;
  size_t SM_nz;
  int *SMp;
  int *SMi;
  double *SMx;
  size_t SM_size = 0;
  if (C == NULL && B->nrow == 0) {
   SM_size = A->ncol;
   SMp = new int[SM_size + 1];
   SMp[0] = 0;
   for (int i = 1; i < A->ncol + 1; ++i) {
    SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]);
   }
  } else if (C == NULL && B->nrow > 0) {
   SM_size = A->ncol + B->nrow;
   SMp = new int[SM_size + 1];
   SMp[0] = 0;
   for (int i = 1; i < A->ncol + 1; ++i) {
    SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
             (B->p[i] - B->p[i - 1]);
   }
  } else if (C->nrow > 0 && B->nrow == 0) {
   SM_size = A->ncol + C->nrow;
   SMp = new int[SM_size + 1];
   SMp[0] = 0;
   for (int i = 1; i < A->ncol + 1; ++i) {
    SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
             (C->p[i] - C->p[i - 1]);
   }
  } else {// both are not null
   SM_size = A->ncol + B->nrow + C->nrow;
   SMp = new int[SM_size + 1];
   SMp[0] = 0;
   for (int i = 1; i < A->ncol + 1; ++i) {
    SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
             (B->p[i] - B->p[i - 1]) +
             (C->p[i] - C->p[i - 1]);
   }
  }
  //Adding diagonal for columns with zero values.
  for (int k = A->ncol + 1; k < SM_size + 1; ++k) {
   SMp[k] = SMp[k - 1] + 1;
  }
  SM_nz = SMp[SM_size];
  SMi = new int[SM_nz];
  SMx = new double[SM_nz]();

  int base1 = A->ncol;
  size_t stp = 0;
  for (int j = 0; j < A->ncol; ++j) {
   stp = SMp[j];
   //Adding Hessian
   for (int i = A->p[j]; i < A->p[j + 1]; ++i) {
    SMi[stp] = A->i[i];
    SMx[stp] = A->x[i];
    stp++;
   }
   base1 = A->ncol;
   //Adding equalities
   if (C != NULL) {
    if (C->nrow > 0) {
     for (int i = C->p[j]; i < C->p[j + 1]; ++i) {
      SMi[stp] = base1 + C->i[i];
      SMx[stp] = C->x[i];
      //std::cout<<"Eq: "<< base1 + C->i[i]<<"; "<<SMx[stp]<<"\n";
      stp++;
     }
    }
    base1 = A->ncol + C->nrow;
   }
   //Adding inequalities
   if (B->nrow > 0) {
    for (int i = B->p[j]; i < B->p[j + 1]; ++i) {
     SMi[stp] = base1 + B->i[i];
     //SMx[stp] = B->x[i];
     stp++;
    }
   }
   assert(stp == SMp[j + 1]);
  }
  //Putting a small value in diagonals
  base1 = A->ncol;
  for (int l = SMp[A->ncol], j = 0; l < SM_nz; ++l, ++j) {
   //SMx[l] = 1e-6;
   SMi[l] = base1 + j;
  }

  SM->ncol = SM->nrow = SM_size;
  SM->p = SMp;
  SM->i = SMi;
  SM->x = SMx;
  SM->stype = -1;
  SM->xtype = CHOLMOD_REAL;
  SM->packed = TRUE;
  SM->sorted = TRUE;
  SM->nzmax = SM_nz;
  sm_rhs = new double[SM->ncol]();
  sm_solution = new double[SM->ncol]();
//  CSC *sKKTt = ptranspose(SM,2,NULL,NULL,0,status);
//  print_csc("skkt\n",SM->ncol,SM->p,SM->i,SM->x);
  //print_csc("\nskkt Trans\n",sKKTt->ncol,sKKTt->p,sKKTt->i,sKKTt->x);
  //std::cout<<"\n";
  return status;
 }

 void SolverSettings::find_perturbation(double tol) {
  for (int i = 0; i < AorSM->ncol; ++i) {
   perturbed_value pv;
   if (std::abs(AorSM->x[AorSM->p[i]]) < tol) {
    pv.col_idx = i;
    pv.per_value = reg_diag;
    AorSM->x[AorSM->p[i]] = reg_diag;
    perturbed_diags.push_back(pv);
   }
  }
  // print_csc("skkt\n",AorSM->ncol,AorSM->p,AorSM->i,AorSM->x);
 }

 void SolverSettings::apply_perturbation(double tol) {
  for (int i = 0; i < A->ncol; ++i) {
   AorSM->x[AorSM->p[i]] += tol;
  }
  for (int i = A->ncol; i < AorSM->ncol; ++i) {
   AorSM->x[AorSM->p[i]] -= tol;
  }
 }

 void SolverSettings::add_perturbation(double tol) {
  for (int i = 0; i < A_ord->ncol; ++i) {
   int ord_col = L->IPerm[i];
   if (ord_col < A->ncol) {
    A_ord->x[A_ord->p[ord_col]] += tol;
   } else {
    A_ord->x[A_ord->p[ord_col]] -= tol;
   }
  }
 }

 void SolverSettings::remove_perturbation(double tol) {
  for (int i = 0; i < A_ord->ncol; ++i) {
   int ord_col = L->IPerm[i];
   if (ord_col < A->ncol) {
    A_ord->x[A_ord->p[ord_col]] -= tol;
   } else {
    A_ord->x[A_ord->p[ord_col]] += tol;
   }
  }
 }

 int SolverSettings::symbolic_analysis() {
  psi->start = psi->tic();

  if (solver_mode == 0) {
   AorSM = A;
   base = 0;
   x = new double[AorSM->ncol];
   if (ldl_variant == 6 || ldl_variant == 7) {
    is_super = 0;
    simplicial_alloc = 1;
   }
  } else { // Mode 1
   is_super = 1; //FIXME: do it in analysis part
   simplicial_alloc = 0;//TODO: this is for low-rank updates
   if (ldl_variant == 6 || ldl_variant == 7
       || ldl_update_variant == 6 || ldl_update_variant == 7) {
    is_super = 0;
    simplicial_alloc = 1;
   }
   build_super_matrix();

   for (int i = 0; i < A->ncol; ++i) {
    sm_rhs[i] = rhs[i];
   }
   rhs = sm_rhs; //FIXME rhs should be read only!
   x = sm_solution;
   AorSM = SM;
   if (C == NULL)
    base = A->ncol;
   else
    base = A->ncol + C->nrow;

  }
  extra_cols = new int[AorSM->ncol]();
#ifdef SYM_REMOV
  if(solver_mode == 1){
    for (int j = base; j < AorSM->ncol; ++j) {
     extra_cols[j] = 1;
    }
   }
#endif
//  for (int j = 147; j < AorSM->ncol; ++j) {
//   extra_cols[j] = 1;
//  }
  //print_vec("Extra: ",0,AorSM->ncol,extra_cols);
  //Fill diagonals with perturbed value
  if (regularization == 1)
   find_perturbation(reg_diag);
  else
   apply_perturbation(reg_diag);
  L = symbolic_analysis_lin_solve(1, AorSM, NULL, NULL, n_relax, z_relax,
                                  AorSM->ncol, prune_ptr, prune_set,
                                  n_level, level_ptr, level_set,
                                  n_par, par_ptr, par_set,
                                  n_level_s, level_ptr_s,
                                  n_par_s, par_ptr_s, par_set_s,
                                  cost_param, level_param, final_seq_node,
                                  status, max_sup_wid, max_col, psi->ordering_time,
                                  simplicial_alloc, extra_cols);
  psi->end = psi->toc();
  psi->analysis_time = psi->elapsed_time(psi->start, psi->end);
  if (L == NULL)
   return 0;
  //
  // print_vec("ordering: ",0,AorSM->ncol,L->Perm);
  //print_vec("inv ordering: ",0,AorSM->ncol,L->IPerm);
  //print_vec("supernode: ",0,L->nsuper,L->super);
  valL = new double[L->xsize]();
  d_val = new double[2 * AorSM->ncol]();
  visible_cnt = new int[L->nsuper]();
  allocate_workspace();
  //ws = new double[2*AorSM->ncol]();
  //ws_int = new int[3*AorSM->ncol]();
  //pinv = new int[AorSM->ncol];
  perm_piv = new int[AorSM->ncol];
  marked = new bool[AorSM->ncol]();
  if (simplicial_alloc) {
   L->x_s = new double[L->xsize_s]();
  }
  AT_ord = ptranspose(AorSM, 2, L->Perm, NULL, 0, status);
  A_ord = ptranspose(AT_ord, 2, NULL, NULL, 0, status);

  etree = L->Parent;
  if (solver_mode == 0) {
   atree = L->sParent;
   if (simplicial_alloc) {
    etree_mod = new int[AorSM->ncol];
    for (int k = 0; k < AorSM->ncol; ++k) {
     etree_mod[k] = L->Parent[k];
    }
    //print_vec("after: \n", 0, AorSM->ncol,etree_mod);
    l_pb = L->p_s; //new int [AorSM->ncol];
    l_pe = NULL; //new int[AorSM->ncol];
    l_i = L->i; //new int[L->xsize]();
    l_x = L->x_s; //new double[L->xsize]();
   }
   //print_vec("\nordering: \n", 0, AorSM->ncol,L->Perm);
  } else { // Mode 1
   //allocate for simplicial factor for low rank update
   l_pb = L->p_s; //new int [AorSM->ncol];
   l_pe = NULL; //new int[AorSM->ncol];
   l_i = L->i; //new int[L->xsize]();
   l_x = L->x_s; //new double[L->xsize]();
   atree = new int[L->nsuper];
   visible_sn = new bool[L->nsuper];
   for (int i = 0; i < L->nsuper; ++i) {
    atree[i] = L->sParent[i];
   }
   for (int l = 0; l < AorSM->ncol; ++l) {
    marked[l] = true;
   }
#ifdef SYM_REMOV
   //Compute the second representation of the tree
    child_sn_ptr = new int[L->nsuper+1];
    child_sn_no = new int[L->nsuper];
    num_sn_child = new int[L->nsuper]();
    populateChildren(L->nsuper,atree,child_sn_ptr,child_sn_no,num_sn_child);
    compressed_set_to_vector(L->nsuper,child_sn_ptr,child_sn_no,children_vec);
    // deleting extra rows symbolically, numerics are already set to zero
    //print_vec("before: \n", 0, L->nsuper,atree);
    //print_vec("\nordering: \n", 0, AorSM->ncol,L->Perm);
    //print_vec("\nsupernode bnds: \n", 0, L->nsuper+1,L->super);
    //print_vec("etree before: \n", 0, AorSM->ncol,L->Parent);
    //int hh = getTreeHeightBruteForce(L->nsuper,atree);
    //print_csc("\noriginal \n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
    //std::cout<<" : "<<hh<<"\n";
    for (int j = base; j < AorSM->ncol; ++j) {//hiding extra cols/rows
     delete_node_tree_simple(j);
    }
   // print_vec("after: \n", 0, L->nsuper,atree);
#endif
   for (int m = 0; m < L->nsuper; ++m) {
    if (marked[m]) {
     modified_sns.push_back(m);
    }
    visible_sn[m] = marked[m];
   }
   architecture_related_params();
  }
  return 1;
 }

 void SolverSettings::reset_symbolic_factor() {
  std::fill(sm_rhs, sm_rhs + SM->ncol, 0);
  std::fill(sm_solution, sm_rhs + SM->ncol, 0);
  std::fill(valL, valL + L->xsize, 0);
  std::fill(d_val, d_val + 2 * AorSM->ncol, 0);
  std::fill(visible_cnt, visible_cnt + L->nsuper, 0);

  std::fill(ws, ws + ws_dbl_size, 0);
  std::fill(ws_int, ws_int + ws_int_size, 0);
  std::fill(ws_zeroed, ws_zeroed + (num_thread * AorSM->ncol), 0);

  //TODO reseting other symbolic info

 }

 void SolverSettings::architecture_related_params() {
  if (num_thread > thread_thresh) {
   s_level_ptr = new int[L->nsuper + 1]();
   s_level_set = new int[L->nsuper]();
   s_level_no = getLevelSet(L->nsuper, L->sParent, s_level_ptr, s_level_set);
  } else {
   s_level_no = -1;
  }
 }

 void SolverSettings::allocate_workspace() {
  //size_t ws_int_size;
  //size_t ws_dbl_size;
  // for update func
  //int workspace: 4*super_max + 2*n + 3*supNo
  // double workspace: 2 * super_max*col_max
  size_t upd_int = 4 * (max_sup_wid + 1) + 2 * AorSM->ncol + 3 * L->nsuper;
  size_t upd_dbl = 2 * (max_sup_wid + 1) * (max_col + 1);

  // temporary used
  size_t tmp_int = 3 * AorSM->ncol;
  size_t tmp_dbl = 2 * AorSM->ncol;

  ws_int_size = std::max(upd_int, tmp_int);
  ws_dbl_size = std::max(upd_dbl, tmp_dbl);

  // ws_size for solve_only= 2*A_ord->ncol
  // ws_size for triangular solves = num_thread*A_ord->ncol
  // ws_size for iter_ref = 4*(max_inner_iter+1) +
  // max_inner_iter*(max_inner_iter+1) + n + max_inner_iter*(n-1) +
  // solve phase
  size_t slve_dbl = 2 * AorSM->ncol;
  if (max_iter > 0) {
   slve_dbl += 4 * (max_inner_iter + 1) +
               (max_inner_iter + 1) * (max_inner_iter + 2) + 1
               + AorSM->ncol +
               (max_inner_iter + 1) * (AorSM->ncol) + 1;
  }
  ws_dbl_size = std::max(slve_dbl, ws_dbl_size);


  //allocating
  ws = new double[ws_dbl_size]();
  ws_int = new int[ws_int_size]();
  ws_zeroed = new double[num_thread * AorSM->ncol]();


 }

 int SolverSettings::numerical_factorization() {
  int ret_val = 0;
  //print_csc("\nORdered: ",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
  switch (ldl_variant) {
   case 1:
//    MKL_Domain_Set_Num_Threads(num_thread, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(num_thread);
    psi->start = psi->tic();
    ret_val = ldl_left_sn_01(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                             L->p, L->s, L->i_ptr, valL,
                             d_val,
                             L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                             atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                             max_sup_wid + 1, max_col + 1, num_pivot);
    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    //MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(1);
    break;
   case 2:
    //MKL_Domain_Set_Num_Threads(num_thread, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(num_thread);
    psi->start = psi->tic();
    ret_val = ldl_left_sn_02_v2(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                                L->p, L->s, L->i_ptr, valL,
                                d_val,
                                L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                max_sup_wid + 1, max_col + 1, num_pivot, perm_piv,
                                L->sParent);
    reorder_matrix();
    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    //MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(1);
    break;
   case 3://parallel static
    psi->start = psi->tic();
#ifdef OPENMP
    ret_val = ldl_left_sn_parallel_01(A_ord->ncol, A_ord->p, A_ord->i,
                                      A_ord->x, L->p, L->s, L->i_ptr, valL,
                                      d_val,
                                      L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                      atree, AT_ord->p, AT_ord->i,
                                      L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                      n_level, level_ptr, level_set,
                                      n_par, par_ptr, par_set,
                                      chunk, num_thread,
                                      max_sup_wid + 1, max_col + 1, num_pivot,
                                      reg_diag);
#endif

    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    break;
   case 4://Parallel SBK
    psi->start = psi->tic();
#ifdef OPENMP
    ret_val = ldl_left_sn_parallel_02(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                                      L->p, L->s, L->i_ptr, valL,
                                      d_val,
                                      L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                      atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                      n_level, level_ptr, level_set,
                                      n_par, par_ptr, par_set,
                                      chunk, num_thread,
                                      max_sup_wid + 1, max_col + 1, num_pivot,
                                      perm_piv);
#endif
    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    reorder_matrix();
    //print_vec("simpl-o: ",0,A_ord->ncol,d_val);
    break;

   case 5://Parallel mixed static SBK
    psi->start = psi->tic();
#ifdef OPENMP
    ret_val = ldl_left_sn_parallel_03(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                                      L->p, L->s, L->i_ptr, valL,
                                      d_val,
                                      L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                      atree, AT_ord->p, AT_ord->i,
                                      L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                      n_level, level_ptr, level_set,
                                      n_par, par_ptr, par_set,
                                      chunk, num_thread,
                                      max_sup_wid + 1, max_col + 1, num_pivot,
                                      perm_piv);
#endif
    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    reorder_matrix();
    break;
   case 6: //Simplicial LDL
    psi->start = psi->tic();
    if (is_super) {
     //convert_supernode_to_simplicial();
     bcsc2csc_aggressive_int(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
                             L->super, valL, l_pb, l_i, l_x);
/*     for (int i = 0; i < L->nzmax; ++i) {
      valL[i]=0;
     }*/
     for (int i = 0; i < A_ord->ncol; ++i) {
      d_val[i] = 0;
     }
     is_super = 0;
    }
    ldl_left_simplicial_02(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                           AT_ord->p, AT_ord->i,
                           l_pb, l_i, l_x, d_val, etree_mod, ws,
                           ws_int);
    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    //print_csc("L:\n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
    break;
   case 7:// parallel simplicial LDL
    psi->start = psi->tic();
    if (is_super) {
     //convert_supernode_to_simplicial();
     bcsc2csc_aggressive_int(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
                             L->super, valL, l_pb, l_i, l_x);
/*     for (int i = 0; i < L->nzmax; ++i) {
      valL[i]=0;
     }*/
     for (int i = 0; i < A_ord->ncol; ++i) {
      d_val[i] = 0;
     }
     is_super = 0;
    }

#ifdef OPENMP
    ldl_parallel_left_simplicial_01(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                                    AT_ord->p, AT_ord->i,
                                    l_pb, l_i, l_x, d_val, etree_mod,
                                    n_level, level_ptr, n_par_s, par_ptr_s,
                                    par_set_s);
#endif
    psi->end = psi->toc();
    psi->fact_time += psi->elapsed_time(psi->start, psi->end);
    break;
   default:
    std::cout << " Wrong algorithm type! \n";
    return -1;
  }
  return 1;
 }

 int SolverSettings::delete_node_tree_simple(int k) {
  // for every child of k, update parent
  int ord_nod = L->IPerm[k];
  int cur_sn = L->col2Sup[ord_nod];
  int cur_par_sn = atree[cur_sn];
  assert(cur_par_sn != INVISIBLE);
  // if it is a root
  if (cur_par_sn == -1) {
   for (int i = 0; i < children_vec[cur_sn].size(); ++i) {
    int cur_c = children_vec[cur_sn][i];
    // only visible nodes
    if (atree[cur_c] != INVISIBLE) {
     //add it to cur_par_sn node
     atree[cur_c] = cur_par_sn;
     //detached_nodes.push_back(cur_c);
    }
   }
   children_vec[cur_sn].clear();
  } else {
   //remove cur_sn from the list of cur_par_sn
   children_vec[cur_par_sn].erase(std::remove(children_vec[cur_par_sn].begin(),
                                              children_vec[cur_par_sn].end(),
                                              cur_sn),
                                  children_vec[cur_par_sn].end());
   assert(cur_sn >= -1);
   for (int i = 0; i < children_vec[cur_sn].size(); ++i) {
    int cur_c = children_vec[cur_sn][i];
    // only visible nodes
    if (atree[cur_c] != INVISIBLE) { // FIXME: no longer needed
     atree[cur_c] = cur_par_sn;
     //add it to cur_par_sn node
     children_vec[cur_par_sn].push_back(cur_c);
    }
   }
  }
  //invalidate supernode correspond to k if all are invisible
  children_vec[cur_sn].clear();
  atree[cur_sn] = INVISIBLE;
  visible_sn[cur_sn] = false;
  marked[cur_sn] = false;
  extra_cols[ord_nod] = 1;
  return 1;
 }

 int SolverSettings::add_node_tree_simple(int k) {
  std::vector<int> tmp;
  // First adjust the parent of k
  int ord_nod = L->IPerm[k];
  int cur_sn = L->col2Sup[ord_nod];
  assert(atree[cur_sn] == INVISIBLE);
  int k_par = L->sParent[cur_sn];
  //up travers to first visible node in original atree
  while (k_par != -1) {
   if (atree[k_par] != INVISIBLE)
    break;
   k_par = L->sParent[k_par];
  }
  assert(children_vec[cur_sn].size() == 0);
  // if the only visible node is a root node
  // add cur_sn to list of k_par children
  if (k_par != ROOT) {
   //children of k_par that has Least common ancestor = cur_sn
   for (int i = 0; i < children_vec[k_par].size(); ++i) {
    int orig_cc = children_vec[k_par][i];
    int cc = children_vec[k_par][i];
    assert(cc != cur_sn);
    while (cc != -1) {
     if (cc == cur_sn) {//cur_sn is LCM
      children_vec[cur_sn].push_back(orig_cc);
      atree[orig_cc] = cur_sn;
      break;
     }
     cc = L->sParent[cc];
    }
    //the node does not have LCM of cur_sn
    if (cc == -1)
     tmp.push_back(orig_cc);
   }
   if (tmp.size() != children_vec[k_par].size()) { // copy the remaining children
    children_vec[k_par].clear();
    if (tmp.size() > 0) {
     children_vec[k_par].insert(children_vec[k_par].begin(), tmp.begin(), tmp.end());
    }
   }
   //add cur_sn to the list of children
   children_vec[k_par].push_back(cur_sn);
  } else {//see if there is any node that cur_sn is the ancestor
   bfs_tree(cur_sn);

/*   for (int i = cur_sn; i >=0 ; --i) {
    for (int j = child_sn_ptr[i]; j < child_sn_ptr[i + 1]; ++j) {
     int cc = child_sn_no[j];
     if(visible_sn[cc] && atree[cc] == ROOT){
      atree[cc] = cur_sn;
      children_vec[cur_sn].push_back(cc);
      if(cc == 18)
       printf("18");
     }
    }
   }*/
  }

  // See if there is any detached node missing here
  // detached nodes can not be an extra col/row
/*  for (int j = 0; j < detached_nodes.size(); ++j) {
   int cc = detached_nodes[j];
   int orig_cc = cc;
   if(atree[cc] == ROOT){ // If it is still a detached node
    assert(cc != cur_sn);
    while(cc != -1){
     if(cc == cur_sn){//cur_sn is LCM
      children_vec[cur_sn].push_back(orig_cc);
      if(orig_cc == 18)
       printf("18");
      atree[orig_cc] = cur_sn;
      break;
     }
     cc = L->sParent[cc];
     if(atree[cc] != INVISIBLE && cc != -1)
      printf("ff");
    }
   }
  }*/
  // modifying the node
  atree[cur_sn] = k_par;
  visible_sn[cur_sn] = true;
  extra_cols[ord_nod] = 0;
  // marked are set outside
  marked[cur_sn] = true;
  // for every original child of k or every node that has k
  // update their parent with cur_sn
  for (int i = child_sn_ptr[cur_sn]; i < child_sn_ptr[cur_sn + 1]; ++i) {
   int cur_c = child_sn_no[i];
   // only visible nodes
   if (atree[cur_c] != INVISIBLE && atree[cur_c] != cur_sn) {
    int par_cur_c = atree[cur_c];
    //update the list of chilren in par_cur_c
    if (par_cur_c != -1) {
     children_vec[par_cur_c].erase(std::remove(children_vec[par_cur_c].begin(),
                                               children_vec[par_cur_c].end(),
                                               cur_c),
                                   children_vec[par_cur_c].end());
    }
    atree[cur_c] = cur_sn;
    children_vec[cur_sn].push_back(cur_c);
   }
  }
  return 1;
 }

 void SolverSettings::bfs_tree(int cur_sn) {
  std::queue<int> qu;
  qu.push(cur_sn);
  while (qu.size() > 0) {
   int vis_node = qu.front();
   qu.pop();
   for (int i = child_sn_ptr[vis_node]; i < child_sn_ptr[vis_node + 1]; ++i) {
    int cc = child_sn_no[i];
    qu.push(cc);
    if (visible_sn[cc] && atree[cc] == ROOT) {
     atree[cc] = cur_sn;
     children_vec[cur_sn].push_back(cc);
    }
   }
  }
 }

 int SolverSettings::delete_node_tree(int k) {
  // for every child of k, update parent
  int ord_nod = L->IPerm[k];
  int par_k = etree_mod[ord_nod];
  int cur_par_sn = L->col2Sup[par_k];
  assert(par_k >= -1);
  for (int i = child_ptr[ord_nod]; i < child_ptr[ord_nod + 1]; ++i) {
   int cur_c = child_no[i];
   // only visible nodes
   if (etree_mod[cur_c] != INVISIBLE) {
    etree_mod[cur_c] = par_k;
    int cur_c_sn = L->col2Sup[cur_c];
    if (cur_c_sn != cur_par_sn) { //Not in the same supernode
     //updating atree
     atree[cur_c_sn] = cur_par_sn;
    }
   }
  }
  etree_mod[ord_nod] = INVISIBLE; //hide the node
  //invalidate supernode correspond to k if all are invisible
  int cur_sn = L->col2Sup[ord_nod];
  visible_cnt[cur_sn]++;
  if (L->super[cur_sn + 1] - L->super[cur_sn] == visible_cnt[cur_sn]) {
   atree[cur_sn] = INVISIBLE;
   marked[cur_sn] = false;
  }
  return 1;
 }

 int SolverSettings::add_node_tree(int k) {
  // First adjust the parent of k
  int ord_nod = L->IPerm[k];
  assert(etree_mod[ord_nod] == INVISIBLE);
  int k_par = etree[ord_nod];
  //up travers to first visible node in original etree
  while (k_par != -1) {
   if (etree_mod[k_par] != INVISIBLE)
    break;
   k_par = etree[k_par];
  }
  etree_mod[ord_nod] = k_par;
  int cur_sn = L->col2Sup[ord_nod];
  // for every child of k, update their parent with k
  for (int i = child_ptr[ord_nod]; i < child_ptr[ord_nod + 1]; ++i) {
   int cur_c = child_no[i];
   // only visible nodes
   if (etree_mod[cur_c] != INVISIBLE) {
    etree_mod[cur_c] = ord_nod;
    int cur_c_sn = L->col2Sup[cur_c];
    if (cur_c_sn != cur_sn) { //Not in the same supernode
     //updating atree
     atree[cur_c_sn] = cur_sn;
    }
   }
  }
  return 1;
 }

 int SolverSettings::update_somod(std::vector<int> add_drp_const, const
 int n_rhs, const double *new_rhs) {
  assert(n_rhs==1); // multiple rhs is not supported yet
  int dim = A_ord->nrow;
  for (int i = 0; i < add_drp_const.size(); ++i) {
   int k = add_drp_const[i];
   for (int j = 0; j < n_rhs; ++j) {
    rhs[(k+base) + j*dim] = new_rhs[k + j*dim];
   }
  }
  return add_del_matrix_qp(1, add_drp_const);
 }

 int SolverSettings::downdate_somod(std::vector<int> add_drp_const, int n_rhs) {
  assert(n_rhs==1); // multiple rhs is not supported yet
  int dim = A_ord->nrow;
  for (int i = 0; i < add_drp_const.size(); ++i) {
   int k = add_drp_const[i];
   for (int j = 0; j < n_rhs; ++j) {
    rhs[(k+base) + j*dim] = 0;
   }
  }
  return add_del_matrix_qp(0, add_drp_const);
 }

 int SolverSettings::add_del_matrix_qp(int add_del, std::vector<int> add_drp_const) {
  int ret_val = 0;
  int *tree;
  int is_simplicial = 0;
  if (simplicial_alloc) {
   is_simplicial = 1;
   tree = etree_mod;
  } else {
   is_simplicial = 0;
   tree = atree;
  }
#ifdef SYM_REMOV
  //Updating symbolic part
   for (int l = 0; l < add_drp_const.size(); ++l) {
    int col = add_drp_const[l];
    //print_vec("before: \n", 0, AorSM->ncol,etree_mod);
    if (add_del) {
     add_node_tree_simple(col+base);
 //    for (int j = BT->p[col]; j < BT->p[col + 1]; ++j) {
 //     int ord_row = BT->i[j];// this is col in kkt matrix
 //     int ord_nod = L->IPerm[ord_row];
 //     int cur_sn = L->col2Sup[ord_nod];
 //     if(atree[cur_sn] == INVISIBLE)
 //      add_node_tree_simple(ord_row);
 //    }
    to_del=0;
    }else{//FIXME
     //delete_node_tree_simple(col+base);
     col_del = col+base;
     to_del=1;
     //printf("RRemoving %d",col);
    }
    //print_vec("after: \n", 0, AorSM->ncol,etree_mod);
   }
#endif
  //print_vec("after: \n", 0, L->nsuper,atree);
  // Numeric part
  /*print_csc("SKKT before adding const\n",A_ord->ncol,A_ord->p,A_ord->i,
            A_ord->x);*/
  for (int l = 0; l < add_drp_const.size(); ++l) {
   int col = add_drp_const[l];
   if (BT->p[col + 1] - BT->p[col] > 0) {
    int ord_col = L->IPerm[base + col]; // this is row in kkt matrix
    /*if(add_del)
     //rhs[base + col] = rhs[base + col];
     printf("TODO");
    else
     rhs[base + col] = 0;*/
    //Finding affected columns
    int s_ord_col = is_simplicial ? ord_col : L->col2Sup[ord_col];
    // For diagonal do it separately, it is not in constraint matrix
    while (s_ord_col >= 0) {
     modified_sns.push_back(s_ord_col);
     s_ord_col = tree[s_ord_col];
    }
/*     if(s_ord_col == -2)
      printf("HRER");*/
    for (int j = BT->p[col]; j < BT->p[col + 1]; ++j) {
     int ord_row = L->IPerm[BT->i[j]];// this is col in kkt matrix
     //traversing etree
     int row_lower, col_lower;
     row_lower = ord_row > ord_col ? ord_row : ord_col;
     col_lower = ord_row > ord_col ? ord_col : ord_row;
     int s_ord_row = is_simplicial ? ord_row : L->col2Sup[ord_row];

     while (s_ord_row >= 0) {
      modified_sns.push_back(s_ord_row);
      s_ord_row = tree[s_ord_row];
     }
/*      if(s_ord_row == -2)
       printf("HRER");*/
     //Search the row in current col of kkt
     for (int k = A_ord->p[col_lower]; k < A_ord->p[col_lower + 1]; ++k) {
      if (A_ord->i[k] == row_lower) {
       if (!add_del) {
        A_ord->x[k] = 0;
       } else {
        assert(A_ord->x[k] == 0);
        A_ord->x[k] = BT->x[j];
       }
       break; // go for the next row of current col.
      }
     }
    }
   }
  }
/*  for (int i = 0; i < L->nsuper; ++i) {
   if(visible_sn[i] &&add_drp_const[0] == 101 && i>1808 && i<1810)
    //if(visible_sn[i] )
    modified_sns.push_back(i);
  }*/
  //print_vec("\n etree :", 0, A_ord->ncol, tree);
//  if(add_drp_const[0] == 47){
//   modified_sns.push_back(156);
//   modified_sns.push_back(158);
//   modified_sns.push_back(159);
//  }
  make_unique(modified_sns);
  for (int k1 = 0; k1 < modified_sns.size(); ++k1) {
   marked[modified_sns[k1]] = 1;
  }
/*    for (int i = 0; i < L->nsuper; ++i) {
   if(visible_sn[i] &&add_drp_const[0] == 219 && i>=0 && i<300)
    atree[i] = L->sParent[i];
  }*/

  //std::cout<<" +: " <<modified_sns.size()<<"\n";

//    if(true){
//   CSC *TMP = ptranspose(A_ord,2,L->IPerm,NULL,0,status);
//   CSC *TMP2 = ptranspose(TMP,2,NULL,NULL,0,status);
//   print_csc("\noriginal \n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
//  }
//  print_vec("PERM: ", 0, A_ord->ncol, L->Perm);
  /*print_csc("\nSKKT after adding const\n",A_ord->ncol,A_ord->p,A_ord->i,
            A_ord->x);*/
  return ret_val;
 }

 int SolverSettings::edit_matrix() {
  return 1;
 }

 int SolverSettings::edit_rhs() {
  return 1;
 }

 double *SolverSettings::edit_solve_rhs() {
  return NULL;
 }

 int SolverSettings::update_factorization() {
  int retval = 0;
  //print_csc("\nORdered: ",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
  switch (ldl_update_variant) {
   case 1:
    //MKL_Domain_Set_Num_Threads(num_thread, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(num_thread);
    psi->start = psi->toc();
    retval = update_ldl_left_sn_01(A_ord->nrow, A_ord->p, A_ord->i, A_ord->x,
                                   L->p, L->s, L->i_ptr, valL,
                                   d_val,
                                   L->super, L->nsuper, psi->timing_chol,
                                   atree, AT_ord->p, AT_ord->i, L->col2Sup,
                                   modified_sns, max_sup_wid + 1, max_col + 1,
                                   num_pivot);
    a_consistent = 1;
    psi->end = psi->toc();
    psi->update_time += psi->elapsed_time(psi->start, psi->end);
    //MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(1);
    break;
   case 2:
    //MKL_Domain_Set_Num_Threads(num_thread, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(num_thread);
    psi->start = psi->toc();
    retval = update_ldl_left_sn_02_v2(A_ord->nrow, A_ord->p, A_ord->i, A_ord->x,
                                      L->p, L->s, L->i_ptr, valL,
                                      d_val,
                                      L->super, L->nsuper, psi->timing_chol,
                                      atree, AT_ord->p, AT_ord->i, L->col2Sup,
                                      modified_sns, max_sup_wid + 1, max_col + 1,
                                      num_pivot, perm_piv, L->sParent, ws_int, ws);
    a_consistent = 0;
    psi->end = psi->toc();
    psi->update_time += psi->elapsed_time(psi->start, psi->end);
    psi->start = psi->toc();
    if (num_pivot > 0)
     reorder_matrix();
    psi->end = psi->toc();
    psi->piv_reord += psi->elapsed_time(psi->start, psi->end);
    //MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
    SET_BLAS_THREAD(1);
    break;
   case 3://parallel static
    std::cout << "Not supported!\n";
    return -1;
    break;
   case 4://Parallel SBK
    psi->start = psi->toc();

#ifdef OPENMP
    retval = update_ldl_left_sn_parallel_02(A_ord->nrow, A_ord->p, A_ord->i, A_ord->x,
                                            L->p, L->s, L->i_ptr, valL,
                                            d_val,
                                            L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                            atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                            n_level, level_ptr, level_set,
                                            n_par, par_ptr, par_set,
                                            chunk, num_thread,
                                            max_sup_wid + 1, max_col + 1, num_pivot,
                                            perm_piv, marked);
#endif

    a_consistent = 0;
    psi->end = psi->toc();
    psi->update_time += psi->elapsed_time(psi->start, psi->end);
    psi->start = psi->toc();
    if (num_pivot > 0)
     reorder_matrix();
    psi->end = psi->toc();
    psi->piv_reord += psi->elapsed_time(psi->start, psi->end);
    break;

   case 5://Parallel mixed static SBK
    std::cout << "Not supported!\n";
    return -1;
    break;
   case 6:
    if (is_super) {
     //convert_supernode_to_simplicial();
     /*bcsc2csc_aggressive_int(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
                             L->super, valL, l_pb, l_i, l_x);*/
     is_super = 0;
    }
    //print_csc("l before: \n",A_ord->ncol,l_pb,l_i,l_x);
    psi->start = psi->toc();

#ifdef OPENMP
    retval = update_ldl_left_simplicial_01(A_ord->ncol, A_ord->p, A_ord->i,
                                           A_ord->x,
                                           AT_ord->p, AT_ord->i,
                                           l_pb, l_i, l_x, d_val,
                                           etree_mod, modified_sns, ws,
                                           ws_int);
#endif
    psi->end = psi->toc();
    psi->update_time += psi->elapsed_time(psi->start, psi->end);
    //print_vec("simpl: ",0,A_ord->ncol,d_val);
    //print_csc("l: \n",A_ord->ncol,l_pb,l_i,l_x);
    break;
   case 7:
    if (is_super) {
     //convert_supernode_to_simplicial();
     /*bcsc2csc_aggressive_int(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
                             L->super, valL, l_pb, l_i, l_x);*/
     is_super = 0;
    }
    //print_vec("mmm : ",0, A_ord->ncol,marked);
    //print_csc("l before: \n",A_ord->ncol,l_pb,l_i,l_x);
    psi->start = psi->toc();

#ifdef OPENMP
    retval = update_ldl_parallel_left_simplicial_01(A_ord->ncol, A_ord->p,
                                                    A_ord->i, A_ord->x,
                                                    AT_ord->p, AT_ord->i,
                                                    l_pb, l_i, l_x, d_val,
                                                    etree_mod, marked,
                                                    n_level, level_ptr, n_par_s,
                                                    par_ptr_s,
                                                    par_set_s);
#endif
    psi->end = psi->toc();
    psi->update_time += psi->elapsed_time(psi->start, psi->end);
    //print_csc("l after: \n",A_ord->ncol,l_pb,l_i,l_x);
    break;
   default:
    std::cout << "Wrong algorithm type!\n";
    return -1;
  }
#ifdef SYM_REMOV
  if(to_del){
    delete_node_tree_simple(col_del);
     }
#endif
  return retval;
 }

 int SolverSettings::edit_update_factorization(int add_del, std::vector<int> add_drp_const) {
  return 1;
 }

 double *SolverSettings::update_solve() {
  return NULL;
 }

 double *SolverSettings::solve_only(const int n_rhs, double *rhs_in) {
  std::copy(rhs_in, rhs_in+A_ord->ncol, rhs);
  return solve_only();
 }

 double *SolverSettings::solve_only() {
  // workspace needed for solve:
  // ws_size for solve_only= 2*A_ord->ncol
  // ws_size for triangular solves = num_thread*A_ord->ncol
  // ws_size for iter_ref = 4*(max_inner_iter+1) +
  // max_inner_iter*(max_inner_iter+1) + n + max_inner_iter*(n-1) +
  //
  psi->start = psi->toc();
  double *x_ord = ws;
  double *rhs_ord = ws + A_ord->ncol;
  for (int i = 0; i < A_ord->ncol; ++i) {
   x_ord[i] = rhs[L->Perm[i]];
  }
  //print_csc("\nA reo: \n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
  //print_vec<double >("\n\nrrhhssss\n: ",0,A_ord->ncol,rhs);
  //print_vec<double >("x reordered: ",0,A->ncol,x_ord);
  //print_vec<int >("ordering: ",0,A_ord->ncol,L->Perm);
  //print_vec<int >("inverse ordering: ",0,A->ncol,L->IPerm);

  if (!is_super) {
   solve_phase_simplicial_ldl(A_ord->ncol, l_pb, l_i, l_x, d_val, x_ord);
  } else if (solver_mode == 1) {
   //std::fill_n(ws+2*A_ord->ncol,num_thread*A_ord->ncol,0);
   solve_phase_ldl_blocked_parallel_permuted_update(A_ord->ncol, d_val, x_ord,
                                                    L->col2Sup, L->super,
                                                    L->p, L->s, valL, L->i_ptr,
                                                    L->nsuper, L->nzmax,
                                                    n_level, level_ptr, level_set,
                                                    n_par, par_ptr, par_set, chunk,
                                                    s_level_no, s_level_ptr, s_level_set,
                                                    visible_sn,
                                                    extra_cols, ws_zeroed);
   //reset all modified cols, it is already updated.
   for (int k1 = 0; k1 < modified_sns.size(); ++k1) {
    marked[modified_sns[k1]] = 0;
   }
   modified_sns.clear();
   //print_vec<double >("x reordered after solve: ",0,A->ncol,x_ord);
   //print_vec<double >("diagonal: ",0,A->ncol,d_val);
  } else {
   solve_phase_ldl_blocked_parallel_permuted(A_ord->ncol, d_val, x_ord,
                                             L->col2Sup, L->super,
                                             L->p, L->s, valL, L->i_ptr,
                                             L->nsuper, L->nzmax,
                                             n_level, level_ptr, level_set,
                                             n_par, par_ptr, par_set, chunk,
                                             L->Perm, L->IPerm);
  }
  if (req_ref_iter > 0) {
   max_iter = req_ref_iter;
   //double *rhs_ord = new double[A_ord->ncol];
   for (int i = 0; i < A_ord->ncol; ++i) {
    rhs_ord[i] = rhs[L->Perm[i]];
   }
   if (regularization == 1)
    set_diags(0);
   else
    remove_perturbation(reg_diag);
   if (solver_mode == 1) {
    //double *ws_nonzero = ws+(2+num_thread)*A_ord->ncol;
    //std::fill_n(ws+2*A_ord->ncol,num_thread*A_ord->ncol,0.0);
    //ws+2*A_ord->ncol is zeroed before
    num_ref_iter = pmgmres_ldlt_auto_update(A_ord->ncol, A_ord->nzmax, A_ord->p, A_ord->i,
                                            A_ord->x,
                                            L->p, L->s, valL, L->nzmax, L->i_ptr,
                                            L->col2Sup, L->super, L->nsuper, d_val,
                                            x_ord, rhs_ord,
                                            max_iter, max_inner_iter,
                                            tol_abs, tol_rel, visible_sn, extra_cols, 1, 2,
                                            n_level, level_ptr, level_set,
                                            n_par, par_ptr, par_set, s_level_no, s_level_ptr,
                                            s_level_set, chunk, ws + 2 * A_ord->ncol,
                                            ws_zeroed);
   } else {
    num_ref_iter = pmgmres_ldlt_auto(A_ord->ncol, A_ord->nzmax, A_ord->p, A_ord->i,
                                     A_ord->x,
                                     L->p, L->s, valL, L->nzmax, L->i_ptr,
                                     L->col2Sup, L->super, L->nsuper, d_val,
                                     x_ord, rhs_ord,
                                     max_iter, max_inner_iter,
                                     tol_abs, tol_rel, 1, 2,
                                     n_level, level_ptr, level_set,
                                     n_par, par_ptr, par_set, chunk);

   }
   //delete []rhs_ord;
   if (regularization == 1)
    set_diags(reg_diag);
   else
    add_perturbation(reg_diag);
  }
  for (int i = 0; i < A_ord->ncol; ++i) {
   x[i] = x_ord[L->IPerm[i]];
  }
  //print_vec<double >("x orig order: ",0,A->ncol,x);
  //print_vec("dval \n",0,2*AorSM->ncol,d_val);
  psi->end = psi->toc();
  psi->solve_time += psi->elapsed_time(psi->start, psi->end);
  return x;
 }

 double *SolverSettings::iterative_ref_only() {
  //double *rhs_ord = new double[A_ord->ncol];
  double *x_ord = ws;
  double *rhs_ord = ws + A_ord->ncol;
  for (int i = 0; i < A_ord->ncol; ++i) {
   double x_tmp = rhs[L->Perm[i]];
   x_ord[i] = x_tmp;
   rhs_ord[i] = x_tmp;
  }
  num_ref_iter = pmgmres_ldlt_cr(A_ord->ncol, A_ord->nzmax, A_ord->p,
                                 A_ord->i,
                                 A_ord->x,
                                 L->p, L->s, valL, L->nzmax, L->i_ptr,
                                 L->col2Sup, L->super, L->nsuper, d_val,
                                 x_ord, rhs_ord,
                                 max_iter, max_inner_iter,
                                 tol_abs, tol_rel, 1, 2,
                                 n_level, level_ptr, level_set,
                                 n_par, par_ptr, par_set, chunk);

  for (int i = 0; i < A_ord->ncol; ++i) {
   x[i] = x_ord[L->IPerm[i]];
  }
  return x;
 }

 void SolverSettings::reorder_matrix() { //TODO: developing a lightweight version of this.

  int *new_ord = ws_int;
  allocateAC(AT_ord, 0, 0, 0, FALSE);
  //print_vec("perm_piv: ",0 ,A_ord->ncol,perm_piv);
  AT_ord = ptranspose(A_ord, 2, perm_piv, NULL, 0, status);
  allocateAC(A_ord, 0, 0, 0, FALSE);
  A_ord = ptranspose(AT_ord, 2, NULL, NULL, 0, status);
  combine_perms(AT_ord->ncol, L->Perm, perm_piv, new_ord);
  for (int i = 0; i < A_ord->ncol; ++i) {
   L->Perm[i] = new_ord[i];
  }
  compute_inv_perm(AT_ord->ncol, L->Perm, L->IPerm);
//  for (int j = 0; j < A_ord->ncol; ++j) {
//   if(perm_piv[j] != j){
//    std::cout<<j<<" -> "<<perm_piv[j]<<";"
//    <<extra_cols[j]<<":"<<extra_cols[perm_piv[j]]<<"\n";
//   }
//  }
//  check_row_idx_l(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
//                      L->super);
//  for (int j = 0; j < A_ord->ncol; ++j) {
//   ws_int[j] = extra_cols[perm_piv[j]];
//  }
//  for (int k = 0; k < A_ord->ncol; ++k) {
//   extra_cols[k] = ws_int[k];
//  }
  //print_csc("\nORdered: ",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
 }

 void SolverSettings::compute_norms() {
  double alp[2] = {1.0, 0};
  double bet[2] = {0.0, 0};
  int norm_type = 0;
  set_diags(0);
  CSC *TMP = ptranspose(A_ord, 2, L->IPerm, NULL, 0, status);
  set_diags(reg_diag);
  CSC *TMP2 = ptranspose(TMP, 2, NULL, NULL, 0, status);
  double *res = new double[TMP2->ncol];
  x_l1 = norm_dense(1, TMP2->ncol, x, norm_type);
  rhs_l1 = norm_dense(1, TMP2->ncol, rhs, norm_type);
  spmv_csc_sym_one_int(TMP2->ncol, TMP2->p, TMP2->i, TMP2->x, -1, alp, bet,
                       1, x, res);
  //print_vec("res mult: ",0,TMP2->ncol,res);
  for (int i = 0; i < TMP2->ncol; ++i) {
   res[i] = rhs[i] - res[i];
  }
  //print_vec("res: ",0,TMP2->ncol,res);
  res_l1 = norm_dense(1, TMP2->ncol, res, norm_type);
  A_l1 = norm_sparse_int(TMP2->ncol, TMP2->p, TMP2->i, TMP2->x, -1, norm_type);
  delete[]res;
  allocateAC(TMP, 0, 0, 0, FALSE);
  allocateAC(TMP2, 0, 0, 0, FALSE);
 }

 double SolverSettings::backward_error() {
  compute_norms();
  if (A_l1 * x_l1 + rhs_l1 > 0)
   bwd_err = res_l1 / (A_l1 * x_l1 + rhs_l1);
  else
   bwd_err = 0;
  std::cout << "d: " << res_l1 << "\n";
  std::cout << "A l1: " << A_l1 << "\n";
  std::cout << "x l1: " << x_l1 << "\n";
  std::cout << "rhs l1: " << rhs_l1 << "\n";
  return bwd_err;
 }

 void SolverSettings::convert_supernode_to_simplicial() {
  size_t actualNNZ = 0;
  for (int i = 0; i < L->nsuper; ++i) {
   int curCol = L->super[i];
   int nxtCol = L->super[i + 1];
   int supWdt = nxtCol - curCol;
   assert(supWdt > 0);
   for (int j = curCol; j < nxtCol; ++j) {
    l_pb[j] = L->p[j] + (j - curCol);
    l_pe[j] = L->p[j + 1];
    for (int k = l_pb[j],
           kk = L->i_ptr[curCol] + (j - curCol);
         k < l_pe[j]; ++k, ++kk) { // copy row indices
     l_i[k] = L->i[kk];
    }
   }
  }
 }

 void SolverSettings::copy_lsuper_to_l() {

 }

 void SolverSettings::set_diags(double d) {
  for (int i = 0; i < perturbed_diags.size(); ++i) {
   int cc = perturbed_diags[i].col_idx;
   int ord_col = L->IPerm[cc];
   double pv = d == 0 ? 0 : perturbed_diags[i].per_value;
   A_ord->x[A_ord->p[ord_col]] = pv;
  }
/*  for (int i = A->ncol; i < A_ord->ncol; ++i) {
   int ord_col = L->IPerm[i];
   //std::cout<<"Diag: "<<A_ord->x[A_ord->p[ord_col]]<<","<<A_ord->i[A_ord->p[ord_col]]<<"\n";
   A_ord->x[A_ord->p[ord_col]] = d;
  }*/
 }

 int SolverSettings::check_ldlt_factor() { // TODO move it to unit test later

  size_t *ia = new size_t[A_ord->ncol + 1]();
  int *ja = new int[L->xsize];
  double *a = new double[L->xsize];
  size_t *A2p = new size_t[A_ord->ncol + 1];
  for (int m = 0; m <= A_ord->ncol; ++m) {
   A2p[m] = static_cast<size_t >(A_ord->p[m]);
  }
  //Converting BCSC to CSC
  bcsc2csc_aggressive(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
                      L->super, valL, ia, ja, a);
  int *tmp = new int[A_ord->ncol + 1]();
  for (int m = 0; m <= A_ord->ncol; ++m) {
   tmp[m] = ia[m];
  }
  print_csc("\n%%MatrixMarket matrix coordinate real general\n",
            A_ord->ncol, tmp, ja, a);
  //print_vec("\nD: ",0,2*A_ord->ncol,d_val);
  print_csc("\n%%MatrixMarket matrix coordinate real symmetric\n",
            A_ord->ncol, A_ord->p, A_ord->i, A_ord->x);
  delete[]tmp;
  double nrm = norm_sparse(A_ord->ncol, ia, ja, a, -1, 1);
  bool check = true; //ldlt_check(A_ord->ncol, ia, ja, a, d_val, A2p,
//    A_ord->i,A_ord->x);
  delete[]ia;
  delete[]ja;
  delete[]a;
  delete[]A2p;

  return check;
 }
}