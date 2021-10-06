//
// Created by kazem on 4/25/19.
//

#ifndef PARS_QP_FORMAT_CONVERTER_H
#define PARS_QP_FORMAT_CONVERTER_H

#include <iostream>
#include <algorithm>
#include "smp_format.h"

namespace format {

/*
 * if a < b return 1
 */
 inline int is_smaller(double a, double b,
                       double tol = std::numeric_limits<double>::min()) {
  if (a - b < tol)
   return 1;
  return 0;
 }

 struct constraint {
  int idx_no;
  double coef; //multiply coef with this number
  constraint() : idx_no(-1), coef(.0) {}
 };

 struct IEForm{
  Description desc;
  CSC *A, *C;
  CSC *AT, *CT;
  Dense *b, *d, *q;
  Dense *duals, *primals;
  CSC *H, *H_general;
  double fixed;
  double optimal_obj;
  IEForm():b(NULLPNTR),d(NULLPNTR),q(NULLPNTR),A(NULLPNTR),AT(NULLPNTR),
           C(NULLPNTR),CT(NULLPNTR),H(NULLPNTR),fixed(0),optimal_obj(0),
           duals(NULLPNTR),primals(NULLPNTR),H_general(NULLPNTR){}

  IEForm(const IEForm *ief){
   A = sym_lib::copy_sparse(ief->A);
   C = sym_lib::copy_sparse(ief->C);
   AT = sym_lib::copy_sparse(ief->AT);
   CT = sym_lib::copy_sparse(ief->CT);
   H = sym_lib::copy_sparse(ief->H);
   H_general = sym_lib::copy_sparse(ief->H_general);
   b = sym_lib::copy_dense(ief->b);
   d = sym_lib::copy_dense(ief->d);
   q = sym_lib::copy_dense(ief->q);
   primals = sym_lib::copy_dense(ief->primals);
   duals = sym_lib::copy_dense(ief->duals);
   optimal_obj = ief->optimal_obj;
   fixed = ief->fixed;
  }

  ~IEForm(){
   delete A;
   delete C;
   delete AT;
   delete CT;
   delete b;
   delete d;
   delete q;
   delete H;
   delete H_general;
   delete duals;
   delete primals;
  }

  bool equality_check(const IEForm* ief, bool is_out = false){
   auto infinity_two_vector = [](Dense *v1, Dense *v2){
    if(!v1 && v2)
     return !v2->is_finite();
    else if(v1 && !v2)
     return !v1->is_finite();
    return sym_lib::are_equal(v1, v2);
   };

   bool h_c = sym_lib::are_equal(H, ief->H);
   auto a_c = sym_lib::are_equal(A, ief->A);
   auto c_c = sym_lib::are_equal(C, ief->C);
   auto b_c = infinity_two_vector(b, ief->b);
   auto d_c = infinity_two_vector(d, ief->d);
   auto q_c = infinity_two_vector(q, ief->q);
   bool pr_c = true, du_c = true, ob_c = true;
   if(is_out){
    pr_c = sym_lib::are_equal(primals, ief->primals);
    du_c = sym_lib::are_equal(duals, ief->duals);
    ob_c = is_equal(optimal_obj, ief->optimal_obj);
   }
   return h_c && a_c && c_c  && d_c && q_c &&
          b_c && b_c && pr_c && du_c && ob_c;
  }


  int get_num_var(){ return H ? H->m : 0;}
  int get_num_eqc(){ return A ? A->m : 0;}
  int get_num_ineqc(){ return C ? C->m : 0;}

  void print(){
   print_csc(H->m, H->n, H->p, H->i, H->x);
   print_csc(A->m, A->n, A->p, A->i, A->x);

   print_csc(C->m, C->n, C->p, C->i, C->x);
   print_dense(d->row,d->col,d->lda,d->a);
  }

 };


 IEForm *load_ie(std::string quad_name, std::string linear_name,
                           std::string eq_name, std::string eql_name,
                           std::string ineq_name, std::string ineql_name){
  auto *ie = new IEForm;
  std::ifstream hin(quad_name);
  if(hin.is_open()){
   read_mtx_csc_real(hin, ie->H);
  }
  hin.close();

  std::ifstream lin(linear_name);
  if(lin.is_open()){
   read_mtx_array_real(lin, ie->q);
  }
  lin.close();

  std::ifstream blin(eql_name);
  if(blin.is_open()){
   read_mtx_array_real(blin, ie->b);
  }
  blin.close();

  std::ifstream Ain(eq_name);
  if(Ain.is_open()){
   read_mtx_csc_real(Ain, ie->A);
  }
  Ain.close();

  std::ifstream Cin(ineq_name);
  if(Cin.is_open()){
   read_mtx_csc_real(Cin, ie->C);
  }
  Cin.close();

  std::ifstream ulin(ineql_name);
  if(ulin.is_open()){
   read_mtx_array_real(ulin, ie->d);
  }
  ulin.close();
 ie->H_general =  sym_lib::make_full(ie->H);

  return ie;
 }



 struct BoundedForm{
  Description desc;
  Dense *l;
  Dense *u;
  Dense *q;
  CSC *A, *AT;
  CSC *H, *H_general;
  Dense *duals, *primals;
  double fixed;
  double optimal_obj;
  BoundedForm():l(NULLPNTR),u(NULLPNTR),q(NULLPNTR),A(NULLPNTR),AT(NULLPNTR),
                H(NULLPNTR),fixed(0),duals(NULLPNTR),primals(NULLPNTR),
                optimal_obj(0), H_general(NULLPNTR){}

  BoundedForm(const BoundedForm *bf){
   A = sym_lib::copy_sparse(bf->A);
   AT = sym_lib::copy_sparse(bf->AT);
   H = sym_lib::copy_sparse(bf->H);
   H_general = sym_lib::copy_sparse(bf->H_general);
   l = sym_lib::copy_dense(bf->l);
   u = sym_lib::copy_dense(bf->u);
   q = sym_lib::copy_dense(bf->q);
   fixed = bf->fixed;
   primals = sym_lib::copy_dense(bf->primals);
   duals = sym_lib::copy_dense(bf->duals);
   optimal_obj = bf->optimal_obj;
  }

  ~BoundedForm(){
   delete A;
   delete AT;
   delete l;
   delete u;
   delete q;
   delete H;
   delete H_general;
   delete duals;
   delete primals;
  }

  bool equality_check(const BoundedForm* bf, bool is_out = false){
   auto infinity_two_vector = [](Dense *v1, Dense *v2){
    if(!v1 && v2)
     return !v2->is_finite();
    else if(v1 && !v2)
     return !v1->is_finite();
    return sym_lib::are_equal(v1, v2);
   };

   bool h_c = sym_lib::are_equal(H, bf->H);
   auto a_c = sym_lib::are_equal(A, bf->A);
   auto l_c = infinity_two_vector(l, bf->l);
   auto u_c = infinity_two_vector(u, bf->u);
   auto q_c = infinity_two_vector(q, bf->q);
   bool p_c = true, d_c = true, o_c = true;
   if(is_out){
    p_c = sym_lib::are_equal(primals, bf->primals);
    d_c = sym_lib::are_equal(duals, bf->duals);
    o_c = is_equal(optimal_obj, bf->optimal_obj);
   }
   return h_c && a_c  && l_c && u_c && q_c &&
          p_c && d_c && o_c;
  }

 };

 BoundedForm *load_bounded(std::string quad_name, std::string linear_name,
   std::string l_name, std::string constraint_name, std::string u_name){
  auto *bf = new BoundedForm;
  std::ifstream hin(quad_name);
  if(hin.is_open()){
   read_mtx_csc_real(hin, bf->H);
  }
  hin.close();

  std::ifstream lin(linear_name);
  if(lin.is_open()){
   read_mtx_array_real(lin, bf->q);
  }
  lin.close();

  std::ifstream blin(l_name);
  if(blin.is_open()){
   read_mtx_array_real(blin, bf->l);
  }
  blin.close();

  std::ifstream Ain(constraint_name);
  if(Ain.is_open()){
   read_mtx_csc_real(Ain, bf->A);
  }
  Ain.close();

  std::ifstream ulin(u_name);
  if(ulin.is_open()){
   read_mtx_array_real(ulin, bf->u);
  }
  ulin.close();
  bf->H_general =  sym_lib::make_full(bf->H);
  return bf;
 }


// for converting smp to ie
 bool find_inequalities_by_bounds(Dense *ld, Dense *ud, CSC *A, CSC *AT,
                                  IEForm* ie_out){
  if(!A){
   ie_out->b = ie_out->d = NULLPNTR;
   ie_out->A = ie_out->AT = ie_out->C = ie_out->CT = NULLPNTR;
   return true;
  }
  bool A_is_transposed = false;
  if(!AT){
   AT = sym_lib::transpose_general(A);
   A_is_transposed = true;
  }
  int num_eq = 0, nnz_eq=0;
  int num_ineq = 0, nnz_ineq=0;
  auto *AT_eq = new CSC;
  auto *AT_ineq = new CSC;
  std::vector<int> eq_idx;
  std::vector<constraint *> ineq_dx;
  int h_dim = A->n;
  int n_const = A->m;
  int *col_cnt_A_eq = new int[h_dim]();
  int *col_cnt_A_ineq = new int[h_dim]();
  double *a_ineq;
  double *l = ld ? ld->a : new double[n_const];
  double *u = ud ? ud->a : new double[n_const];
  if(!ld)
   std::fill_n(l, n_const, MIN_DBL);
  if(!ud)
   std::fill_n(l, n_const, MAX_DBL);
  //ie_out->b = new Dense(2*n_const,1,1);; a_eq = ie_out->b->a;
  ie_out->d = new Dense(2*n_const,1,1);; a_ineq = ie_out->d->a;
  ie_out->AT = AT_eq;
  ie_out->CT = AT_ineq;
  for (int i = 0; i < A->m; ++i) {
   l[i] = std::abs(l[i]) < 1e-14 ? 0 : l[i];
   u[i] = std::abs(u[i]) < 1e-14 ? 0 : u[i];

   // Invalid constraint
   if ((is_equal(l[i], MIN_DBL) && is_equal(u[i], MAX_DBL)) ||
       (is_equal(l[i], MAX_DBL) && is_equal(u[i], MAX_DBL)) ||
       (is_equal(l[i], MIN_DBL) && is_equal(u[i], MIN_DBL))) {
    continue;
   }
   if (is_equal(l[i], u[i])) {//eq
    num_eq++;
    nnz_eq += (AT->p[i + 1] - AT->p[i]);
    eq_idx.push_back(i);
    for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
     col_cnt_A_eq[AT->i[j]]++;
    }
   } else { // ineq
    constraint *c_csnt = new constraint;
    c_csnt->idx_no = i;
    if (is_equal(l[i], MIN_DBL)) {//one constraint Ax<=b
     a_ineq[num_ineq] = u[i];
     num_ineq++;
     nnz_ineq += (AT->p[i + 1] - AT->p[i]);
     c_csnt->coef = 1;
     ineq_dx.push_back(c_csnt);
     for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
      col_cnt_A_ineq[AT->i[j]]++;
     }
    } else if (is_equal(u[i], MAX_DBL)) {//one constraint Ax>=b ==> Ax <= -b
     a_ineq[num_ineq] = -l[i];
     num_ineq++;
     nnz_ineq += (AT->p[i + 1] - AT->p[i]);
     c_csnt->coef = -1;
     ineq_dx.push_back(c_csnt);
     for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
      col_cnt_A_ineq[AT->i[j]]++;
     }
    } else {//two constraints
     a_ineq[num_ineq] = u[i];
     a_ineq[num_ineq + 1] = -l[i];
     num_ineq += 2;
     nnz_ineq += 2 * (AT->p[i + 1] - AT->p[i]);
     c_csnt->coef = 1;
     ineq_dx.push_back(c_csnt);
     constraint *c_csnt_2 = new constraint;
     c_csnt_2->idx_no = i;
     c_csnt_2->coef = -1;
     ineq_dx.push_back(c_csnt_2);
     for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
      col_cnt_A_ineq[AT->i[j]] += 2;
     }
    }

   }
  }
  assert(num_ineq <= 2 * A->m);
  assert(num_eq == 0); // all equalities assumed to be sperated
  /// building the equality constraint matrix
  /*if(num_eq > 0) {
   AT_eq->m = h_dim;
   AT_eq->n = num_eq;
   AT_eq->p = new int[num_eq + 1];
   AT_eq->i = new int[nnz_eq];
   AT_eq->x = new double[nnz_eq];
   assert(eq_idx.size() == num_eq);
   AT_eq->p[0] = 0;
   for (int ll = 0; ll < eq_idx.size(); ++ll) {
    int cur_idx = eq_idx[ll];
    int nnz_cur_row = AT->p[cur_idx + 1] - AT->p[cur_idx];
    AT_eq->p[ll + 1] = AT_eq->p[ll] + nnz_cur_row;
    for (int ii = AT_eq->p[ll], jj = AT->p[cur_idx];
         ii < AT_eq->p[ll + 1]; ++ii, ++jj) {
     AT_eq->i[ii] = AT->i[jj];
     AT_eq->x[ii] = AT->x[jj];
    }
   }
   ie_out->b->row = num_eq;
   ie_out->A = sym_lib::transpose_general(AT_eq);
  }else{
   delete ie_out->b;
   ie_out->b = NULLPNTR;
   delete AT_eq;
   ie_out->A = ie_out->AT = NULLPNTR;
  }*/

  //building the ineq constraint matrix
  if(num_ineq>0){
   AT_ineq->m = h_dim;
   AT_ineq->n = num_ineq;
   AT_ineq->p = new int[num_ineq + 1];
   AT_ineq->i = new int[nnz_ineq];
   AT_ineq->x = new double[nnz_ineq];
   assert(ineq_dx.size() == num_ineq);
   AT_ineq->p[0] = 0;
   for (int ll = 0; ll < ineq_dx.size(); ++ll) {
    int cur_idx = ineq_dx[ll]->idx_no;
    double cur_coef = ineq_dx[ll]->coef;
    int nnz_cur_row = AT->p[cur_idx + 1] - AT->p[cur_idx];
    AT_ineq->p[ll + 1] = AT_ineq->p[ll] + nnz_cur_row;
    for (int ii = AT_ineq->p[ll], jj = AT->p[cur_idx];
         ii < AT_ineq->p[ll + 1]; ++ii, ++jj) {
     AT_ineq->i[ii] = AT->i[jj];
     AT_ineq->x[ii] = cur_coef * AT->x[jj];
    }
   }
   AT_ineq->nnz = AT_ineq->p[num_ineq];
   ie_out->C = sym_lib::transpose_general(AT_ineq);
   ie_out->d->row = num_ineq;
  } else{
   delete ie_out->d;
   ie_out->d = NULLPNTR;
   delete AT_ineq;
   ie_out->C = ie_out->CT = NULLPNTR;
  }

  eq_idx.clear();
  for (int ii = 0; ii < ineq_dx.size(); ++ii) {
   delete ineq_dx[ii];
  }
  if(!ld)
   delete []l;
  if(!ud)
   delete []u;
  delete[]col_cnt_A_ineq;
  delete[]col_cnt_A_eq;
  if(A_is_transposed)
   delete AT;
  return true;
 }

 /// Used for converting bounded to SMP
 /// \param ld
 /// \param ud
 /// \param A
 /// \param AT
 /// \param smp_out
 /// \return
 bool find_smp_inequalities_by_bounds(Dense *ld, Dense *ud, CSC *A,
                                      CSC *AT, SMP* smp_out){
  if(!A){
   smp_out->b_ = smp_out->l_ = smp_out->u_ = NULLPNTR;
   smp_out->A_ = smp_out->AT_ = smp_out->C_ = smp_out->CT_ = NULLPNTR;
   return true;
  }
  bool A_is_transposed = false;
  if(!AT){
   AT = sym_lib::transpose_general(A);
   A_is_transposed = true;
  }
  int num_eq = 0, nnz_eq=0;
  int num_ineq = 0, nnz_ineq=0;
  auto *AT_eq = new CSC;
  auto *AT_ineq = new CSC;
  std::vector<int> eq_idx;
  std::vector<constraint *> ineq_dx;
  int h_dim = A->n;
  int n_const = A->m;
  int *col_cnt_A_eq = new int[h_dim]();
  int *col_cnt_A_ineq = new int[h_dim]();
  double *a_eq;
  double *l_ineq;
  double *u_ineq;
  double *l = ld ? ld->a : new double[n_const];
  double *u = ud ? ud->a : new double[n_const];
  if(!ld)
   std::fill_n(l, n_const, MIN_DBL);
  if(!ud)
   std::fill_n(l, n_const, MAX_DBL);
  smp_out->b_ = new Dense(n_const,1,1); a_eq = smp_out->b_->a;
  smp_out->l_ = new Dense(n_const,1,1); l_ineq = smp_out->l_->a;
  smp_out->u_ = new Dense(n_const,1,1); u_ineq = smp_out->u_->a;
  smp_out->CT_ = AT_ineq;
  for (int i = 0; i < A->m; ++i) {
   if ((is_equal(l[i], MIN_DBL) && is_equal(u[i], MAX_DBL)) ||
       (is_equal(l[i], MAX_DBL) && is_equal(u[i], MAX_DBL)) ||
       (is_equal(l[i], MIN_DBL) && is_equal(u[i], MIN_DBL))) {
    continue;
   }
   if (is_equal(l[i], u[i])) {//eq
    a_eq[num_eq] = l[i];
    num_eq++;
    nnz_eq += (AT->p[i + 1] - AT->p[i]);
    eq_idx.push_back(i);
    for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
     col_cnt_A_eq[AT->i[j]]++;
    }
   } else { // ineq
    auto *c_csnt = new constraint;
    c_csnt->idx_no = i;
    l_ineq[num_ineq] = l[i];
    u_ineq[num_ineq] = u[i];
    num_ineq++;
    nnz_ineq += (AT->p[i + 1] - AT->p[i]);
    c_csnt->coef = 1;
    ineq_dx.push_back(c_csnt);
    for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
     col_cnt_A_ineq[AT->i[j]]++;
    }
   }
  }
  assert(num_ineq <=  A->m);
  assert(num_eq <= A->m);
  /// building the equality constraint matrix
  if(num_eq > 0){
   AT_eq->m = h_dim;
   AT_eq->n = num_eq;
   AT_eq->p = new int[num_eq + 1];
   AT_eq->i = new int[nnz_eq];
   AT_eq->x = new double[nnz_eq];

   assert(eq_idx.size() == num_eq);
   AT_eq->p[0] = 0;
   for (int ll = 0; ll < eq_idx.size(); ++ll) {
    int cur_idx = eq_idx[ll];
    int nnz_cur_row = AT->p[cur_idx + 1] - AT->p[cur_idx];
    AT_eq->p[ll + 1] = AT_eq->p[ll] + nnz_cur_row;
    for (int ii = AT_eq->p[ll], jj = AT->p[cur_idx];
         ii < AT_eq->p[ll + 1]; ++ii, ++jj) {
     AT_eq->i[ii] = AT->i[jj];
     AT_eq->x[ii] = AT->x[jj];
    }
   }
   smp_out->b_->row = num_eq;
   AT_eq->nnz = nnz_eq;
   assert(AT_eq->p[num_eq] == nnz_eq);
  smp_out->AT_ = AT_eq;
  smp_out->A_ = sym_lib::transpose_general(AT_eq);
  } else{
   delete smp_out->b_;
   delete AT_eq;
   smp_out->AT_ = smp_out->A_ = NULLPNTR;
   smp_out->b_ = NULLPNTR;
  }

  ///building the ineq constraint matrix
  if(num_ineq >0){
   AT_ineq->m = h_dim;
   AT_ineq->n = num_ineq;
   AT_ineq->p = new int[num_ineq + 1];
   AT_ineq->i = new int[nnz_ineq];
   AT_ineq->x = new double[nnz_ineq];
   assert(ineq_dx.size() == num_ineq);
   AT_ineq->p[0] = 0;
   for (int ll = 0; ll < ineq_dx.size(); ++ll) {
    int cur_idx = ineq_dx[ll]->idx_no;
    double cur_coef = ineq_dx[ll]->coef;
    int nnz_cur_row = AT->p[cur_idx + 1] - AT->p[cur_idx];
    AT_ineq->p[ll + 1] = AT_ineq->p[ll] + nnz_cur_row;
    for (int ii = AT_ineq->p[ll], jj = AT->p[cur_idx];
         ii < AT_ineq->p[ll + 1]; ++ii, ++jj) {
     AT_ineq->i[ii] = AT->i[jj];
     AT_ineq->x[ii] = cur_coef * AT->x[jj];
    }
   }
   AT_ineq->nnz = AT_ineq->p[num_ineq];
   smp_out->C_ = sym_lib::transpose_general(AT_ineq);
   smp_out->l_->row = smp_out->u_->row = num_ineq;
  } else{
   delete smp_out->l_;
   delete smp_out->u_;
   delete AT_ineq;
   smp_out->CT_ = smp_out->C_ = NULLPNTR;
   smp_out->l_ = smp_out->u_ = NULLPNTR;
  }

  eq_idx.clear();
  for (int ii = 0; ii < ineq_dx.size(); ++ii) {
   delete ineq_dx[ii];
  }
  if(!ld)
   delete []l;
  if(!ud)
   delete []u;
  delete[]col_cnt_A_ineq;
  delete[]col_cnt_A_eq;
  if(A_is_transposed)
   delete AT;
  return true;
 }



 struct QPFormatConverter {
  std::string problem_name;
  std::string desc_;

  // Bounded format l <= A <= u
  bool bounded_converted;
  BoundedForm *bf_;

  double *l;
  double *u;
  CSC *A, *AT;
  CSC *H_full;

  // Inequality/Equality (IE) format
  bool ie_converted;
  IEForm *ief_;
  CSC *A_eq, *A_ineq;
  CSC *AT_eq, *AT_ineq;
  double *a_eq, *a_ineq;
  CSC *H;

  // for QL testing
  bool dense_converted;
  int ql_wanted;
  double *AB_d, *ab_eqineq, *H_d;
  double *A_d, *B_d;

  // SMP format
  bool smp_converted;
  format::SMP *smp_;

  // linear term
  double *q; // common in all formats
  int mode; //0: IE format given, 1: bounded given, 2: ? 3: ?
  double min_dbl, max_dbl;

  int num_eq, num_ineq;
  int nnz_eq, nnz_ineq;

  QPFormatConverter():smp_converted(false), ie_converted(false),
                      bounded_converted(false), dense_converted(false), smp_(NULLPNTR),
                      ief_(NULLPNTR), bf_(NULLPNTR),
                      l(NULLPNTR), u(NULLPNTR), A(NULLPNTR), AT(NULLPNTR), H_full(NULLPNTR),
                      A_eq(NULLPNTR), A_ineq(NULLPNTR), AT_eq(NULLPNTR), AT_ineq(NULLPNTR),
                      a_eq(NULLPNTR), a_ineq(NULLPNTR), H(NULLPNTR), AB_d(NULLPNTR), ab_eqineq(NULLPNTR),
                      H_d(NULLPNTR), A_d(NULLPNTR), B_d(NULLPNTR), q(NULLPNTR){
   mode = 0;
   max_dbl = MAX_DBL;//std::numeric_limits<double >::max();
   min_dbl = -MAX_DBL;//std::numeric_limits<double >::min();
   num_eq = 0;
   num_ineq = 0;
   nnz_eq = 0;
   nnz_ineq = 0;
   problem_name = "noname";
   ql_wanted = 0;
   A_d = B_d = NULL;
  }


  QPFormatConverter(CSC *H_full_in, double *q_in, CSC *A_in, double *l_in, double *u_in) :
    H_full(H_full_in), q(q_in), A(A_in), l(l_in), u(u_in) {
   mode = 1;
   max_dbl = MAX_DBL;//std::numeric_limits<double >::max();
   min_dbl = -MAX_DBL;//std::numeric_limits<double >::min();
   num_eq = 0;
   num_ineq = 0;
   nnz_eq = 0;
   nnz_ineq = 0;
   problem_name = "noname";
   ql_wanted = 0;
   A_d = B_d = NULL;

  }

  explicit QPFormatConverter(BoundedForm *bf):QPFormatConverter(){
   if(bf){
    bf_ = new BoundedForm(bf);
    bounded_converted = true;
   }
  }

  explicit QPFormatConverter(IEForm *ief):QPFormatConverter(){
   if(ief){
    ief_ = new IEForm(ief);
    ie_converted = true;
   }
  }

  explicit QPFormatConverter(SMP *smp):QPFormatConverter(){
   if(smp){
    smp_ = new SMP(smp);
    smp_converted = true;
   }
  }

  ~QPFormatConverter() {
   delete smp_;
   delete bf_;
   delete ief_;
  }

  bool smp_to_ie(){
   if(ie_converted)
    return true;
   if(!smp_converted)
    return false;
   ief_ = new IEForm;
   ief_->desc = smp_->desc_struct_;
   problem_name = ief_->desc.name_;
   ief_->H = sym_lib::copy_sparse(smp_->H_);
   ief_->q = sym_lib::copy_dense(smp_->q_);
   ief_->fixed = smp_->r_;
   ief_->b = sym_lib::copy_dense(smp_->b_);
   ief_->A = sym_lib::copy_sparse(smp_->A_);
   ief_->AT = sym_lib::transpose_general(ief_->A);
   if(!smp_->l_){ //no specific conversion is required, almost there
    ief_->C = sym_lib::copy_sparse(smp_->C_);
    ief_->d = sym_lib::copy_dense(smp_->u_);
    ief_->CT = sym_lib::transpose_general(ief_->C);
   } else{// converting lower bound to upperbound
    find_inequalities_by_bounds(smp_->l_, smp_->u_, smp_->C_,ief_->CT, ief_);
   }
   // To avoid null pointer needed for the benchmark
   if(!ief_->A){
    ief_->A = new CSC(0, num_var(),0, false,GENERAL);
   }
   if(!ief_->C){
    ief_->C = new CSC(0, num_var(),0, false,GENERAL);
   }
   ief_->H_general = sym_lib::make_full(ief_->H);
   ie_converted = true;
   return true;
  }

  bool ie_to_smp(){
   if(smp_converted)// smp is already there.
    return true;
   if(!ie_converted) //ie is neither loaded nor converted
    return false;
   int num_vars = ief_->H->n;
   smp_ = new SMP("");
   smp_->desc_struct_=ief_->desc;
   smp_->desc_=smp_->desc_struct_.get_desc();
   problem_name = smp_->desc_struct_.name_;
   smp_->H_ = sym_lib::copy_sparse(ief_->H);
   smp_->H_->stype = LOWER;
   smp_->q_ = sym_lib::copy_dense(ief_->q);
   smp_->r_ = ief_->fixed;
   smp_->A_ = sym_lib::copy_sparse(ief_->A);
   smp_->A_->stype = GENERAL;
   smp_->b_ = sym_lib::copy_dense(ief_->b);
   smp_->l_ = NULLPNTR;
   smp_->C_ = sym_lib::copy_sparse(ief_->C);
   smp_->C_->stype = GENERAL;
   smp_->u_ = sym_lib::copy_dense(ief_->d);
   smp_converted = true;
   return true;
  }

  bool bounded_to_smp(){
   if(smp_converted)// smp is already there.
    return true;
   if(!bounded_converted) //bounded is neither loaded nor converted
    return false;
   smp_ = new SMP("");
   smp_->desc_struct_=bf_->desc;
   smp_->desc_=smp_->desc_struct_.get_desc();
   problem_name = smp_->desc_struct_.name_;
   smp_->H_ = sym_lib::copy_sparse(bf_->H);
   smp_->H_->stype = LOWER;
   smp_->q_ = sym_lib::copy_dense(bf_->q);
   smp_->r_ = bf_->fixed;
   find_smp_inequalities_by_bounds(bf_->l, bf_->u, bf_->A, bf_->AT, smp_);
   smp_->A_->stype = GENERAL;

   return true;
  }

  bool smp_to_bounded(){
   if(bounded_converted)// bounded is already there.
    return true;
   if(!smp_converted) //smp is neither loaded nor converted
    return false;
   int num_vars = smp_->H_->n;
   int n_eq = smp_->A_ ? smp_->A_->m : 0;
   int n_ineq = smp_->C_ ? smp_->C_->m : 0;
   bf_ = new BoundedForm;
   bf_->desc = smp_->desc_struct_;
   bf_->H = sym_lib::copy_sparse(smp_->H_);
   bf_->H->stype = LOWER;
   bf_->H_general = sym_lib::make_full(bf_->H);
   bf_->q = sym_lib::copy_dense(smp_->q_);
   bf_->fixed = smp_->r_;
   bf_->A = sym_lib::concatenate_two_CSC(smp_->A_, smp_->C_);
   bf_->A->stype=GENERAL;
   if((!smp_->l_ && smp_->u_) || (smp_->l_ && !smp_->u_)){
    if(!smp_->l_){
     auto *d = new Dense(smp_->u_->row, 1, 1, min_dbl);
     bf_->l = sym_lib::concatenate_two_dense(smp_->b_,d);
     bf_->u = sym_lib::concatenate_two_dense(smp_->b_,smp_->u_);
     delete d;
    }
    if(!smp_->u_){
     auto *d = new Dense(smp_->l_->row, 1, 1, max_dbl);
     bf_->l = sym_lib::concatenate_two_dense(smp_->b_,smp_->l_);
     bf_->u = sym_lib::concatenate_two_dense(smp_->b_,d);
     delete d;
    }
   }else{
    bf_->l = sym_lib::concatenate_two_dense(smp_->b_,smp_->l_);
    bf_->u = sym_lib::concatenate_two_dense(smp_->b_,smp_->u_);
   }
   //bf_->l = sym_lib::copy_dense(smp_->l_);
   // To avoid null pointer for benchmark
   if(!bf_->l && bf_->A) bf_->l = new Dense(bf_->A->m, 1, 1, min_dbl);
   if(!bf_->u && bf_->A) bf_->u = new Dense(bf_->A->m, 1, 1, max_dbl);
   if(!bf_->A)
    bf_->A = new CSC(0, num_var(),0, false,GENERAL);

   bounded_converted = true;
   return true;
  }

  bool load_smp(const std::string& smp_file){
   if(smp_converted)
    return false;
   smp_ = new SMP(smp_file);
   if(smp_->load()){
    smp_converted = true;
    return true;
   }
   return false;
  }

  size_t num_ineq_constraints(){
   if(smp_converted)
    return smp_->C_ ? smp_->C_->m : 0;
   if(ie_converted)
    return ief_->C ? ief_->C->m : 0;
   if(bounded_converted){
    bounded_to_smp();
    return smp_->C_ ? smp_->C_->m : 0;
   }
   return 0;
  }

  size_t num_var(){
   if(smp_converted)
    return smp_->H_ ? smp_->H_->m : 0;
   if(ie_converted)
    return ief_->H ? ief_->H->m : 0;
   if(bounded_converted){
    return bf_->H ? bf_->H->m : 0;
   }
   return 0;
  }

  void print_log() {
   if (ief_) {
    std::cout << std::setprecision(20) << ief_->desc.name_ << "," << ief_->get_num_var() << "," << ief_->H->nnz << ",";
    std::cout << std::setprecision(20) << ief_->get_num_eqc() << "," << (ief_->A ? ief_->A->nnz : 0) <<
    "," << ief_->get_num_ineqc() << ",";
    std::cout << std::setprecision(20) << (ief_->C ? ief_->C->nnz : 0) << ",";
   }
  }

  /*int read_IE_format(std::string hessian_file,
                     std::string linear_file,
                     std::string eq_file,
                     std::string eq_file_u,
                     std::string ineq_file,
                     std::string ineq_file_u) {
   problem_name = linear_file;
   H = new CSC;
   CSC *H_tmp = new CSC;
   A_eq = new CSC;
   A_ineq = new CSC;
   AT_eq = new CSC;
   AT_ineq = new CSC;
   if (!readMatrix(hessian_file, H_tmp->ncol, H_tmp->nzmax, H_tmp->p,
                   H_tmp->i, H_tmp->x))
    return 0;
   //print_csc("HH\n ",H_tmp->ncol,H_tmp->p,H_tmp->i,H_tmp->x);
   bool is_expanded = expandMatrix(H_tmp->ncol, H_tmp->nzmax, H_tmp->p,
                                   H_tmp->i, H_tmp->x,
                                   H->nzmax, H->p,
                                   H->i, H->x);
   if (is_expanded) {
    //H has the expanded version already.
    H->ncol = H->nrow = H_tmp->ncol;
    allocateAC(H_tmp, 0, 0, 0, FALSE);
   } else {
    H->x = H_tmp->x;
    H->p = H_tmp->p;
    H->i = H_tmp->i;
    H->nzmax = H_tmp->nzmax;
    H->ncol = H->nrow = H_tmp->ncol;
   }
   H->nzmax = H->p[H->ncol];
   H->stype = -1;
   H->packed = 1;
   //TODO: if it is full symmetric, make it half
*//*  CSC *lower_H = computeLowerTriangular(H);
  if(lower_H){
   allocateAC(H_tmp,0,0,0,FALSE);
   H->ncol = H->nrow = lower_H->ncol;
   H->x = lower_H->x;
   H->p = lower_H->p;
   H->i = lower_H->i;
   H->nzmax = lower_H->nzmax;
  }*//*

   if (eq_file != "none") {
    if (!readMatrix_rect(eq_file, A_eq->nrow, A_eq->ncol, A_eq->nzmax,
                         A_eq->p, A_eq->i, A_eq->x))
     return 0;
    if (A_eq->nrow > 0) {
     A_eq->nzmax = A_eq->p[A_eq->ncol];
     num_eq = A_eq->nrow;
     a_eq = new double[A_eq->nrow];
     read_vector(eq_file_u, A_eq->nrow, a_eq);
     transpose_unsym(A_eq->nrow, A_eq->ncol, A_eq->p, A_eq->i,
                     A_eq->x, AT_eq->nrow, AT_eq->ncol, AT_eq->p, AT_eq->i,
                     AT_eq->x);
    } else {
     A_eq->nzmax = A_eq->nrow = A_eq->ncol = 0;
     A_eq->p = A_eq->i = NULL;
     A_eq->x = a_eq = NULL;
     num_eq = 0;
    }
   } else {
    A_eq->nzmax = A_eq->nrow = A_eq->ncol = 0;
    A_eq->p = A_eq->i = NULL;
    A_eq->x = a_eq = NULL;
    num_eq = 0;
   }

   if (ineq_file != "none") {
    if (!readMatrix_rect(ineq_file, A_ineq->nrow, A_ineq->ncol, A_ineq->nzmax,
                         A_ineq->p, A_ineq->i, A_ineq->x))
     return 0;
    if (A_ineq->nrow > 0) {
     A_ineq->nzmax = A_ineq->p[A_ineq->ncol];
     num_ineq = A_ineq->nrow;

     a_ineq = new double[A_ineq->nrow];
     read_vector(ineq_file_u, A_ineq->nrow, a_ineq);
     transpose_unsym(A_ineq->nrow, A_ineq->ncol, A_ineq->p, A_ineq->i,
                     A_ineq->x, AT_ineq->nrow, AT_ineq->ncol, AT_ineq->p,
                     AT_ineq->i, AT_ineq->x);
    } else {
     A_ineq->nzmax = A_ineq->nrow = A_ineq->ncol = 0;
     A_ineq->p = A_ineq->i = NULL;
     A_ineq->x = a_ineq = NULL;
     num_ineq = 0;
    }
   } else {
    A_ineq->nzmax = A_ineq->nrow = A_ineq->ncol = 0;
    A_ineq->p = A_ineq->i = NULL;
    A_ineq->x = a_ineq = NULL;
    num_ineq = 0;
   }

   q = new double[H->ncol];
   read_vector(linear_file, H->ncol, q);

   mode = 0;
   delete H_tmp;
   return 1;
  }

  int read_bounded_format(std::string hessian_file,
                          std::string linear_file,
                          std::string ineq_file_l,
                          std::string ineq_file,
                          std::string ineq_file_u) {
   problem_name = linear_file;
   H_full = new CSC;
   CSC *H_tmp = new CSC;
   A = new CSC;
   if (!readMatrix(hessian_file, H_tmp->ncol, H_tmp->nzmax, H_tmp->p,
                   H_tmp->i, H_tmp->x))
    return 0;
   //print_csc("HH\n ",H_tmp->ncol,H_tmp->p,H_tmp->i,H_tmp->x);
   bool is_expanded = expandMatrix(H_tmp->ncol, H_tmp->nzmax, H_tmp->p,
                                   H_tmp->i, H_tmp->x,
                                   H_full->nzmax, H_full->p,
                                   H_full->i, H_full->x);
   if (is_expanded) {
    H_full->ncol = H_full->nrow = H_tmp->ncol;
    allocateAC(H_tmp, 0, 0, 0, FALSE);
   } else {
    H_full->x = H_tmp->x;
    H_full->p = H_tmp->p;
    H_full->i = H_tmp->i;
    H_full->nzmax = H_tmp->nzmax;
    H_full->ncol = H_full->nrow = H_tmp->ncol;
    delete H_tmp;
   }

   if (!readMatrix_rect(ineq_file, A->nrow, A->ncol, A->nzmax, A->p, A->i, A->x))
    return 0;
   H_full->nzmax = H_full->p[H_full->ncol];
   A->nzmax = A->p[A->ncol];
   q = new double[H_full->ncol];
   read_vector(linear_file, H_full->ncol, q);
   l = new double[A->nrow];
   read_vector(ineq_file_l, A->nrow, l);
   u = new double[A->nrow];
   read_vector(ineq_file_u, A->nrow, u);
   mode = 3;
   return 1;
  }

  int print_bounded_format() {
   print_csc("\nObjective function: \n", H_full->ncol, H_full->p, H_full->i,
             H_full->x);
   print_csc("\nConstraint matrix: \n", A->ncol, A->p, A->i, A->x);
   print_vec("\nLower bounds: \n", 0, A->nrow, l);
   print_vec("\nUpper bounds: \n", 0, A->nrow, u);
  }

  int B2IE() {
   if (mode != 1 && mode != 3)
    return 0;
   // Objective function
   H = new CSC;
   make_lower(H_full->ncol, H_full->nzmax, H_full->p, H_full->i, H_full->x,
              H->ncol, H->nzmax, H->p, H->i, H->x);

   // Constraints
   AT = new CSC;
   transpose_unsym(A->nrow, A->ncol, A->p, A->i, A->x, AT->nrow, AT->ncol,
                   AT->p, AT->i, AT->x);
   //print_csc("\nAT: \n",AT->ncol,AT->p,AT->i,AT->x);
   A_eq = new CSC;
   A_ineq = new CSC;
   AT_eq = new CSC;
   AT_ineq = new CSC;
   std::vector<int> eq_idx;
   std::vector<constraint *> ineq_dx;
   int h_dim = H_full->ncol;
   int *col_cnt_A_eq = new int[h_dim]();
   int *col_cnt_A_ineq = new int[h_dim]();
   a_eq = new double[2 * A->nrow];//FIXME: allocate it smarter
   a_ineq = new double[2 * A->nrow];

   for (int i = 0; i < A->nrow; ++i) {
*//*   l[i] = std::abs(l[i]) < 1e-14 ? 0 : l[i];
   u[i] = std::abs(u[i]) < 1e-14 ? 0 : u[i];*//*
    // Invalid constraint
    if ((is_equal(l[i], min_dbl) && is_equal(u[i], max_dbl)) ||
        (is_equal(l[i], max_dbl) && is_equal(u[i], max_dbl)) ||
        (is_equal(l[i], min_dbl) && is_equal(u[i], min_dbl))) {
     continue;
    }
    if (is_equal(l[i], u[i])) {//eq

     a_eq[num_eq] = l[i];
     num_eq++;
     nnz_eq += (AT->p[i + 1] - AT->p[i]);
     eq_idx.push_back(i);
     for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
      col_cnt_A_eq[AT->i[j]]++;
     }
    } else { // ineq
     constraint *c_csnt = new constraint;
     c_csnt->idx_no = i;
     if (is_equal(l[i], min_dbl)) {//one constraint Ax<=b
      a_ineq[num_ineq] = u[i];
      num_ineq++;
      nnz_ineq += (AT->p[i + 1] - AT->p[i]);
      c_csnt->coef = 1;
      ineq_dx.push_back(c_csnt);
      for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
       col_cnt_A_ineq[AT->i[j]]++;
      }
     } else if (is_equal(u[i], max_dbl)) {//one constraint Ax>=b ==> Ax <= -b
      a_ineq[num_ineq] = -l[i];
      num_ineq++;
      nnz_ineq += (AT->p[i + 1] - AT->p[i]);
      c_csnt->coef = -1;
      ineq_dx.push_back(c_csnt);
      for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
       col_cnt_A_ineq[AT->i[j]]++;
      }
     } else {//two constraints
      a_ineq[num_ineq] = u[i];
      a_ineq[num_ineq + 1] = -l[i];
      num_ineq += 2;
      nnz_ineq += 2 * (AT->p[i + 1] - AT->p[i]);
      c_csnt->coef = 1;
      ineq_dx.push_back(c_csnt);
      constraint *c_csnt_2 = new constraint;
      c_csnt_2->idx_no = i;
      c_csnt_2->coef = -1;
      ineq_dx.push_back(c_csnt_2);
      for (int j = AT->p[i]; j < AT->p[i + 1]; ++j) {
       col_cnt_A_ineq[AT->i[j]] += 2;
      }
     }

    }
   }
   assert(num_ineq <= 2 * A->nrow);
   assert(num_eq <= 2 * A->nrow);
   // building the equality constraint matrix
   AT_eq->nrow = h_dim;
   AT_eq->ncol = num_eq;
   AT_eq->p = new int[num_eq + 1];
   AT_eq->i = new int[nnz_eq];
   AT_eq->x = new double[nnz_eq];
   assert(eq_idx.size() == num_eq);
   AT_eq->p[0] = 0;
   for (int ll = 0; ll < eq_idx.size(); ++ll) {
    int cur_idx = eq_idx[ll];
    int nnz_cur_row = AT->p[cur_idx + 1] - AT->p[cur_idx];
    AT_eq->p[ll + 1] = AT_eq->p[ll] + nnz_cur_row;
    for (int ii = AT_eq->p[ll], jj = AT->p[cur_idx];
         ii < AT_eq->p[ll + 1]; ++ii, ++jj) {
     AT_eq->i[ii] = AT->i[jj];
     AT_eq->x[ii] = AT->x[jj];
    }
   }
   transpose_unsym(AT_eq->nrow, AT_eq->ncol, AT_eq->p, AT_eq->i,
                   AT_eq->x, A_eq->nrow, A_eq->ncol, A_eq->p, A_eq->i, A_eq->x);
   A_eq->nzmax = nnz_eq;
   //building the ineq constraint matrix
   AT_ineq->nrow = h_dim;
   AT_ineq->ncol = num_ineq;
   AT_ineq->p = new int[num_ineq + 1];
   AT_ineq->i = new int[nnz_ineq];
   AT_ineq->x = new double[nnz_ineq];
   assert(ineq_dx.size() == num_ineq);
   AT_ineq->p[0] = 0;
   for (int ll = 0; ll < ineq_dx.size(); ++ll) {
    int cur_idx = ineq_dx[ll]->idx_no;
    double cur_coef = ineq_dx[ll]->coef;
    int nnz_cur_row = AT->p[cur_idx + 1] - AT->p[cur_idx];
    AT_ineq->p[ll + 1] = AT_ineq->p[ll] + nnz_cur_row;
    for (int ii = AT_ineq->p[ll], jj = AT->p[cur_idx];
         ii < AT_ineq->p[ll + 1]; ++ii, ++jj) {
     AT_ineq->i[ii] = AT->i[jj];
     AT_ineq->x[ii] = cur_coef * AT->x[jj];
    }
   }

   transpose_unsym(AT_ineq->nrow, AT_ineq->ncol, AT_ineq->p, AT_ineq->i,
                   AT_ineq->x, A_ineq->nrow, A_ineq->ncol, A_ineq->p,
                   A_ineq->i, A_ineq->x);
   A_ineq->nzmax = nnz_ineq;
   A_eq->stype = 0;
   A_ineq->stype = 0;
   H_full->stype = 0;
   eq_idx.clear();
   for (int ii = 0; ii < ineq_dx.size(); ++ii) {
    delete ineq_dx[ii];
   }
   delete[]col_cnt_A_ineq;
   delete[]col_cnt_A_eq;
  }

  int print_IE_format() {

   print_csc("\nObjective function: \n", H->ncol, H->p, H->i,
             H->x);
   print_vec("\nq: \n", 0, H->ncol, q);
   print_csc("\nEq matrix: \n", A_eq->ncol, A_eq->p, A_eq->i, A_eq->x);
   print_csc("\nIneq matrix: \n", A_ineq->ncol, A_ineq->p, A_ineq->i, A_ineq->x);
   print_vec("\nEq bounds: \n", 0, num_eq, a_eq);
   print_vec("\nIneq bounds: \n", 0, num_ineq, a_ineq);
  }

  int IE_export_to_dense(std::string mat_name) {

   if (mode != 1 && mode != 3 && mode != 0)
    return 0;

   int total_cnst = A_ineq->nrow + A_eq->nrow;
   if (ql_wanted) {
    if (mode != 0) {
     AB_d = new double[(total_cnst) * H_full->ncol]();
     ab_eqineq = new double[total_cnst];
    } else {
     H_full = new CSC;
     int status = 0;
     CSC *HT = ptranspose(H, 2, NULL, NULL, 0, status);
     int nnzH = H->nzmax, nnz_full = 0;
     make_full(H->ncol, nnzH, H->p, H->i, H->x, HT->p, HT->i, HT->x, nnz_full,
               H_full->p, H_full->i, H_full->x);
     H_full->nrow = H_full->ncol = H->ncol;
     H_full->stype = 0;
     H_full->packed = 1;
     AB_d = new double[(total_cnst) * H_full->ncol]();
     ab_eqineq = new double[total_cnst];
     allocateAC(HT, 0, 0, 0, FALSE);
    }
   }
   H_d = new double[H_full->ncol * H_full->ncol]();
   sparse2dense(H_full, H_d);
   //write_dense(mat_name+"_P.dns",H_full->nrow,H_full->ncol,H_d);
   //write_vector(mat_name+"_q.dns",H_full->ncol,q);
   if (A_eq->nrow > 0) {
    A_d = new double[A_eq->nrow * A_eq->ncol]();
    sparse2dense(A_eq, A_d);
    if (ql_wanted) {
     for (int i = 0; i < A_eq->ncol; ++i) {
      for (int j = 0; j < A_eq->nrow; ++j) {
       AB_d[i * total_cnst + j] = -A_d[i * A_eq->nrow + j];
      }
     }
     for (int k = 0; k < A_eq->nrow; ++k) {
      ab_eqineq[k] = a_eq[k];
     }
    }
    //write_dense(mat_name+"_A.dns",A_eq->nrow,A_eq->ncol,A_d);
    //write_vector(mat_name+"_aeq.dns",A_eq->nrow,a_eq);
    //delete []A_d;
   }
   if (A_ineq->nrow > 0) {
    B_d = new double[A_ineq->nrow * A_ineq->ncol]();
    sparse2dense(A_ineq, B_d);
    if (ql_wanted) {
     for (int i = 0; i < A_ineq->ncol; ++i) {
      for (int j = 0; j < A_ineq->nrow; ++j) {
       AB_d[i * total_cnst + j + A_eq->nrow] = -B_d[i * A_ineq->nrow + j];
      }
     }
     for (int k = 0; k < A_ineq->nrow; ++k) {
      ab_eqineq[k + A_eq->nrow] = a_ineq[k];
     }
    }
    // write_dense(mat_name+"_B.dns",A_ineq->nrow,A_ineq->ncol,B_d);
#if 0
    double *tmp=NULL;
    int n_row_tmp,n_col_tmp;
    read_dense(mat_name+"_B.dns",n_row_tmp,n_col_tmp,tmp);
    if(n_row_tmp!=A_ineq->nrow || n_col_tmp != A_ineq->ncol)
     std::cout<<"wrong dims\n";
    for (int i = 0; i < n_row_tmp * n_col_tmp; ++i) {
     if(tmp[i]-B_d[i]>0){
      std::cout<<"Wrong dense conversion\n";
     }
    }
    delete []tmp;
#endif
    // write_vector(mat_name+"_bineq.dns",A_ineq->nrow,a_ineq);
#if 0
    double *tmp_v = new double[A_ineq->nrow];
    read_vector(mat_name+"_bineq.dns",A_ineq->nrow,tmp_v);
    for (int i = 0; i < A_ineq->nrow; ++i) {
     if(tmp_v[i]-a_ineq[i]>0){
      std::cout<<"Wrong dense conversion\n";
     }
    }
    delete []tmp_v;
#endif
    //delete []B_d;
   }
   if (!ql_wanted)
    delete[]H_d;
  }

  void print_log() {
   std::setprecision(48);
   if (mode == 1 || mode == 3 || mode == 0) {
    std::cout << problem_name << "," << H->ncol << "," << H->nzmax << ",";
    std::cout << A_eq->nrow << "," << A_eq->nzmax << "," << A_ineq->nrow << ",";
    std::cout << A_ineq->nzmax << ",";
   }
  }*/

 };

}
#endif //PARS_QP_FORMAT_CONVERTER_H
