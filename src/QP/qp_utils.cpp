//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/QP/qp_utils.h"

#include "nasoq/common/Norm.h"
#include "nasoq/common/Util.h"
#include "nasoq/matrixVector/spmv_CSC.h"

#include <algorithm>


namespace nasoq {

 double dot1(int n, double *a, double *b) {
  double result = 0.0;
  for (int i = 0; i < n; ++i) {
   result += (a[i] * b[i]);
  }
  return result;
 }

 int
 build_qp_from_file(std::string hessian_file, std::string ineq_file, std::string ineq_file_rhs, std::string linear_file,
                    std::string opt_primal_file, std::string opt_dual_file, std::string opt_obj_file,
                    std::string warm_start_file, size_t &sizeH, size_t &nnzH, double *&q, int *&colH, int *&rowH,
                    double *&valH, size_t &ad1, size_t &ad2, size_t &nnzA, int *&colA, int *&rowA, double *&valA,
                    double *&b_ineq, double *&optimal_primal, double &objective, double *&optimal_dual, double *&init_x) {

  if (!readMatrix(hessian_file, sizeH, nnzH, colH, rowH, valH))
   return 0;
  if (!readMatrix_rect(ineq_file, ad1, ad2, nnzA, colA, rowA, valA))
   return 0;

  q = new double[sizeH];
  read_vector(linear_file, sizeH, q);
  b_ineq = new double[ad1];
  read_vector(ineq_file_rhs, ad1, b_ineq);

  optimal_primal = new double[sizeH];
  read_vector(opt_primal_file, sizeH, optimal_primal);
  optimal_dual = new double[ad1];
  if (opt_dual_file != "none") {
   read_vector(opt_dual_file, ad1, optimal_dual);
  }
  if (opt_obj_file != "none") {
   double *obj_tmp = new double[1];
   read_vector(opt_obj_file, 1, obj_tmp);
   objective = obj_tmp[0];
   delete[]obj_tmp;
  }

  if (warm_start_file != "none") {
   init_x = new double[sizeH];
   read_vector(warm_start_file, sizeH, init_x);
  }
  return 1;
 }

 int build_qp_optimality_from_file(std::string hessian_file, std::string ineq_file, std::string ineq_file_rhs,
                                   std::string linear_file, std::string opt_primal_file, std::string opt_dual_file,
                                   size_t &sizeH, size_t &nnzH, double *&q, int *&colH, int *&rowH, double *&valH,
                                   size_t &ad1, size_t &ad2, size_t &nnzA, int *&colA, int *&rowA, double *&valA,
                                   double *&b_ineq, double *&optimal_primal, double *&optimal_dual) {

  if (!readMatrix(hessian_file, sizeH, nnzH, colH, rowH, valH))
   return 0;
  if (!readMatrix_rect(ineq_file, ad1, ad2, nnzA, colA, rowA, valA))
   return 0;

  q = new double[sizeH];
  read_vector(linear_file, sizeH, q);
  b_ineq = new double[ad1];
  read_vector(ineq_file_rhs, ad1, b_ineq);

  optimal_primal = new double[sizeH];
  read_vector(opt_primal_file, sizeH, optimal_primal);
  optimal_dual = new double[ad1];
  if (opt_dual_file != "none") {
   read_vector(opt_dual_file, ad1, optimal_dual);
  }

  return 1;
 }

 void
 build_qp_01(size_t &sizeH, size_t &nnzH, double *&q, int *&colH, int *&rowH, double *&valH, size_t &ad1, size_t &ad2,
             size_t &nnzA, int *&colA, int *&rowA, double *&valA, double *&b_ineq, double *&optimal_primal,
             double &objective) {

  sizeH = 2;
  nnzH = 3;
  ad1 = 4;
  ad2 = 2;
  nnzA = 8;
  q = new double[sizeH];
  colH = new int[sizeH + 1];
  rowH = new int[nnzH];
  valH = new double[nnzH];
  colA = new int[ad2 + 1];
  rowA = new int[nnzA];
  valA = new double[nnzA];
  b_ineq = new double[ad1];
  optimal_primal = new double[sizeH];

  q[0] = -4;
  q[1] = -4;

  colH[0] = 0;
  colH[1] = 2;
  colH[2] = 4;
  rowH[0] = 0;
  rowH[1] = 1;
  rowH[2] = 1;
  valH[0] = 2;
  valH[1] = 0;
  valH[2] = 2;

  colA[0] = 0;
  colA[1] = 4;
  colA[2] = 8;
  rowA[0] = 0;
  rowA[1] = 1;
  rowA[2] = 2;
  rowA[3] = 3;
  rowA[4] = 0;
  rowA[5] = 1;
  rowA[6] = 2;
  rowA[7] = 3;
  valA[0] = 2;
  valA[1] = 1;
  valA[2] = -1;
  valA[3] = -2;
  valA[4] = 1;
  valA[5] = -1;
  valA[6] = -1;
  valA[7] = 1;
  b_ineq[0] = 2;
  b_ineq[1] = 1;
  b_ineq[2] = 1;
  b_ineq[3] = 2;

  optimal_primal[0] = 0.4;
  optimal_primal[1] = 1.2;
 }

 void
 build_qp_02(size_t &sizeH, size_t &nnzH, double *&q, int *&colH, int *&rowH, double *&valH, size_t &ad1, size_t &ad2,
             size_t &nnzA, int *&colA, int *&rowA, double *&valA, double *&b_ineq, double *&optimal_primal,
             double &objective) {

  sizeH = 2;
  nnzH = 3;
  ad1 = 5;
  ad2 = 2;
  nnzA = 10;
  q = new double[sizeH];
  colH = new int[sizeH + 1];
  rowH = new int[nnzH];
  valH = new double[nnzH];
  colA = new int[ad2 + 1];
  rowA = new int[nnzA];
  valA = new double[nnzA];
  b_ineq = new double[ad1];
  optimal_primal = new double[sizeH];

  q[0] = -2;
  q[1] = -5;

  colH[0] = 0;
  colH[1] = 2;
  colH[2] = 3;
  rowH[0] = 0;
  rowH[1] = 1;
  rowH[2] = 1;
  valH[0] = 2;
  valH[1] = 0;
  valH[2] = 2;

  colA[0] = 0;
  colA[1] = 5;
  colA[2] = 10;
  rowA[0] = 0;
  rowA[1] = 1;
  rowA[2] = 2;
  rowA[3] = 3;
  rowA[4] = 4;
  rowA[5] = 0;
  rowA[6] = 1;
  rowA[7] = 2;
  rowA[8] = 3;
  rowA[9] = 4;

  valA[0] = -1;
  valA[1] = 1;
  valA[2] = 1;
  valA[3] = -1;
  valA[4] = 0;
  valA[5] = 2;
  valA[6] = 2;
  valA[7] = -2;
  valA[8] = 0;
  valA[9] = -1;
  b_ineq[0] = 2;
  b_ineq[1] = 6;
  b_ineq[2] = 2;
  b_ineq[3] = 0;
  b_ineq[4] = 0;

  optimal_primal[0] = 1.4;
  optimal_primal[1] = 1.7;
 }

 void
 build_qp_03(size_t &sizeH, size_t &nnzH, double *&q, int *&colH, int *&rowH, double *&valH, size_t &ad1, size_t &ad2,
             size_t &nnzA, int *&colA, int *&rowA, double *&valA, double *&b_ineq, double *&optimal_primal,
             double &objective) {

  sizeH = 3;
  nnzH = 6; //lower part
  ad1 = 3;
  ad2 = 3;
  nnzA = 9;
  q = new double[sizeH];
  colH = new int[sizeH + 1];
  rowH = new int[nnzH];
  valH = new double[nnzH];
  colA = new int[ad2 + 1];
  rowA = new int[nnzA];
  valA = new double[nnzA];
  b_ineq = new double[ad1];
  optimal_primal = new double[sizeH];

  q[0] = 6;
  q[1] = 0;
  q[2] = 2;
/*
 *
 8, -2, 1,
 -2, 8, 4,
 1, 4, 8
 */
  colH[0] = 0;
  colH[1] = 3;
  colH[2] = 5;
  colH[3] = 6;
  rowH[0] = 0;
  rowH[1] = 1;
  rowH[2] = 2;
  rowH[3] = 1;
  rowH[4] = 2;
  rowH[5] = 2;

  valH[0] = 8;
  valH[1] = -2;
  valH[2] = 1;
  valH[3] = 8;
  valH[4] = 4;
  valH[5] = 8;

  colA[0] = 0;
  colA[1] = 3;
  colA[2] = 6;
  colA[3] = 9;
  rowA[0] = 0;
  rowA[1] = 1;
  rowA[2] = 2;
  rowA[3] = 0;
  rowA[4] = 1;
  rowA[5] = 2;
  rowA[6] = 0;
  rowA[7] = 1;
  rowA[8] = 2;

  valA[0] = -1;
  valA[1] = 0;
  valA[2] = -1;
  valA[3] = 0;
  valA[4] = -1;
  valA[5] = -1;
  valA[6] = -2;
  valA[7] = 0;
  valA[8] = -3;
  b_ineq[0] = 0;
  b_ineq[1] = 0;
  b_ineq[2] = -2;

  //-0.5135135135134224  0.                  0.8378378378377807
  optimal_primal[0] = -0.5135135135135135;
  optimal_primal[1] = 0;
  optimal_primal[2] = 0.8378378378378379;
  objective = 2.02703;
 }

 void
 build_qp_04(size_t &sizeH, size_t &nnzH, double *&q, int *&colH, int *&rowH, double *&valH, size_t &ad1, size_t &ad2,
             size_t &nnzA, int *&colA, int *&rowA, double *&valA, double *&b_ineq, size_t &bd1, size_t &bd2,
             size_t &nnzB, int *&colB, int *&rowB, double *&valB, double *&b_eq, double *&optimal_primal,
             double &objective) {

  sizeH = 3;
  nnzH = 6; //lower part
  ad1 = 3;
  ad2 = 3;
  nnzA = 9;
  bd1 = 1;
  bd2 = 3;
  nnzB = 3;
  q = new double[sizeH];
  colH = new int[sizeH + 1];
  rowH = new int[nnzH];
  valH = new double[nnzH];

  colA = new int[ad2 + 1];
  rowA = new int[nnzA];
  valA = new double[nnzA];
  b_ineq = new double[ad1];

  colB = new int[bd2 + 1];
  rowB = new int[nnzB];
  valB = new double[nnzB];
  b_eq = new double[bd1];

  optimal_primal = new double[sizeH];

  q[0] = 6;
  q[1] = 0;
  q[2] = 2;
/*
 *
 8, -2, 1,
 -2, 8, 4,
 1, 4, 8
 */
  colH[0] = 0;
  colH[1] = 3;
  colH[2] = 5;
  colH[3] = 6;
  rowH[0] = 0;
  rowH[1] = 1;
  rowH[2] = 2;
  rowH[3] = 1;
  rowH[4] = 2;
  rowH[5] = 2;

  valH[0] = 8;
  valH[1] = -2;
  valH[2] = 1;
  valH[3] = 8;
  valH[4] = 4;
  valH[5] = 8;

  colA[0] = 0;
  colA[1] = 3;
  colA[2] = 6;
  colA[3] = 9;
  rowA[0] = 0;
  rowA[1] = 1;
  rowA[2] = 2;
  rowA[3] = 0;
  rowA[4] = 1;
  rowA[5] = 2;
  rowA[6] = 0;
  rowA[7] = 1;
  rowA[8] = 2;

  valA[0] = -1;
  valA[1] = 0;
  valA[2] = -1;
  valA[3] = 0;
  valA[4] = -1;
  valA[5] = -1;
  valA[6] = -2;
  valA[7] = 0;
  valA[8] = -3;
  b_ineq[0] = 0;
  b_ineq[1] = 0;
  b_ineq[2] = -2;

  colB[0] = 0;
  colB[1] = 1;
  colB[2] = 2;
  colB[3] = 3;
  rowB[0] = 0;
  rowB[1] = 0;
  rowB[2] = 0;

  valB[0] = 0;
  valB[1] = 0;
  valB[2] = 2;
  b_eq[0] = 1;

  optimal_primal[0] = 2.50e-02;
  optimal_primal[1] = 4.75e-01;
  optimal_primal[2] = 5.00e-01;
  objective = 3.9937500251437976;
 }

 bool CMP(idx_val a, idx_val b) {
  return a.val < b.val;
 }

 int remove_from_list(std::vector<idx_val> &list_const, int c_no) {
  std::vector<int>::iterator it;
  for (int i = 0; i < list_const.size(); ++i) {
   if (list_const[i].idx == c_no) {
    list_const.erase(list_const.begin() + i);
    return 1;
   }
  }
  return 0;
 }

 double compute_primal_obj(double *primal_vars, CSC *H, double *q) {
  double primal_obj = 0;
  int status = 0;
  double alp[2] = {1.0, 0};
  double bet[2] = {0.0, 0};
  double *tmp = new double[H->ncol];
  spmv_csc_sym_one_int(H->ncol, H->p, H->i, H->x, -1, alp, bet,
                       1, primal_vars, tmp);
  /*CSC *HTT = ptranspose(H,2,NULL,NULL,0,status);
  primal_obj = quad_form(HTT,primal_vars);*/
  primal_obj = 0.5 * dot1(H->ncol, tmp, primal_vars);
  primal_obj += dot1(H->ncol, q, primal_vars);
  delete[]tmp;
  return primal_obj;
 }

 double constraint_sat_norm(CSC *B, CSC *A, double *b, double *a, double *primal_vars) {
  double cons_sat_norm = 0;
  double *tmp = new double[B->nrow + A->nrow]();
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
  for (int i = 0; i < B->nrow; ++i) {
   double diff = tmp[i] - b[i];
   //tmp[i] = std::min(tmp[i]-b[i],0.0);
   //tmp[i] = std::max(b[i] - tmp[i], 0.0); //Equal to top
   tmp[i] = diff < 0.0 ? 0.0 : diff;
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
   //cons_sat_norm += constraint_sat_norm_eq;
   cons_sat_norm = std::max(cons_sat_norm, constraint_sat_norm_eq);
  }
  delete[]tmp;
  return cons_sat_norm;
 }

 double
 lagrangian_residual_norm(CSC *H, CSC *B, CSC *BT, CSC *A, CSC *AT, double *q, double *primal_vars, double *dual_vars,
                          double *dual_vars_eq) {
  double lag_res = 0;
  double alp[2] = {1.0, 0};
  double bet[2] = {0.0, 0};
  double *workspace = new double[4 * H->ncol]();
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

  //print_vec("Beq_dual",0,BT->ncol,dual_vars);
  if (B->nrow > 0) {
   spmv_csc_small(BT->nrow, BT->ncol, BT->p, BT->i, BT->x,
                  dual_vars, Bdual);
  } else {
   for (int i = 0; i < H->ncol; ++i) {
    Bdual[i] = 0;
   }
  }
  //print_vec("Bdual",0,H->ncol,Bdual);
  for (int i = 0; i < H->ncol; ++i) {
   res_norm[i] = Hx[i] + q[i] + Bdual[i] + B_eq_dual[i];
  }
  //print_vec("res_norm",0,H->ncol,res_norm);
  lag_res = norm_dense(H->ncol, 1, res_norm, 0);
  delete[]workspace;
  return lag_res;
 }

 double complementarity_norm(CSC *B, CSC *A, double *b, double *a, double *primal_vars, double *dual_vars) {
  double complementarity_infn = 0;
  double *tmp = new double[B->nrow + A->nrow]();
  spmv_csc_small(B->nrow, B->ncol, B->p, B->i, B->x, primal_vars, tmp);
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] = dual_vars[i] * (tmp[i] - b[i]);
  }
  complementarity_infn = norm_dense(B->nrow, 1, tmp, 0);
  delete[]tmp;
  return complementarity_infn;
 }

 double non_negativity_norm(CSC *B, double *dual_vars) {
  double *tmp = new double[B->nrow];
  for (int i = 0; i < B->nrow; ++i) {
   tmp[i] = dual_vars[i] > 0.0 ? 0 : dual_vars[i];
  }
  double non_negativity_infn = norm_dense(B->nrow, 1, tmp, 0);
  delete []tmp;
  return non_negativity_infn;
 }

 void Fischer_Burmeister_func(int n, double *a, double *b, double *f) {
  for (int i = 0; i < n; ++i) {
   double a_i = a[i];
   double b_i = b[i];
   f[i] = a_i + b_i - sqrt(a_i * a_i + b_i * b_i);
  }
 }

 void bound_to_general_ineq() {

 }

 qp_info::qp_info() : num_refactor(0), num_solve(0), num_update(0), num_downdate(0),
                      factt(0.0), tot(0.0), init(0.0), sw(true) {}

 std::chrono::time_point<std::chrono::system_clock> qp_info::tic() {
  return std::chrono::system_clock::now();
 }

 std::chrono::time_point<std::chrono::system_clock> qp_info::toc() {
  return std::chrono::system_clock::now();
 }

 double qp_info::elapsed_time(std::chrono::time_point<std::chrono::system_clock> beg,
                              std::chrono::time_point<std::chrono::system_clock> lst) {
  double ret = 0;
  elapsed_seconds = lst - beg;
  ret = elapsed_seconds.count();
  return ret;
 }

 void qp_info::print() {
  std::cout << "refact, update, downdate, solve: " << num_refactor <<
            ", " << num_update << ", " << num_downdate << ", " << num_solve << "\n";
 }

 idx_val::idx_val() : val(0), idx(0) {
 }

 idx_val::idx_val(int i, double v) : val(v), idx(i) {
 }
}