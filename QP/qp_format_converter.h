//
// Created by kazem on 4/25/19.
//

#ifndef PARS_QP_FORMAT_CONVERTER_H
#define PARS_QP_FORMAT_CONVERTER_H

#include <iostream>
#include <algorithm>
#include "def.h"
#include "Util.h"
#include "SparseUtils.h"
#include "Transpose.h"



inline int is_equal(double a, double b,
  double tol=std::numeric_limits<double>::min()){
 if(std::abs(a - b) <= tol)
  return 1;
 return 0;
}

/*
 * if a < b return 1
 */
inline int is_smaller(double a, double b,
                    double tol=std::numeric_limits<double>::min()){
 if(a - b < tol)
  return 1;
 return 0;
}

struct constraint{
 int idx_no;
 double coef; //multiply coef with this number
 constraint():idx_no(-1),coef(.0){}
};

struct QPFormatConverter{

 std::string problem_name;
 // Bounded format l <= A <= u
 double *l;
 double *u;
 CSC *A, *AT;
 CSC *H_full;

 // Inequality/Equality (IE) format
 CSC *A_eq,*A_ineq;
 CSC *AT_eq,*AT_ineq;
 double *a_eq, *a_ineq;
 CSC *H;

 // for QL testing
 int ql_wanted;
 double *AB_d, *ab_eqineq, *H_d;
 double  *A_d, *B_d;
 // linear term
 double *q; // common in all formats
 int mode; //0: IE format given, 1: bounded given, 2: ? 3: ?
 double min_dbl, max_dbl;

 int num_eq, num_ineq;
 int nnz_eq, nnz_ineq;

 QPFormatConverter(){
  mode=0;
  max_dbl = 1e+20;//std::numeric_limits<double >::max();
  min_dbl = -1e+20;//std::numeric_limits<double >::min();
  num_eq=0; num_ineq=0;
  nnz_eq=0; nnz_ineq=0;
  problem_name = "noname";
  ql_wanted = 0;
  A_d = B_d = NULL;
 }

 QPFormatConverter(CSC *H_full_in, double *q_in, CSC *A_in, double *l_in, double *u_in):
 H_full(H_full_in),q(q_in), A(A_in),l(l_in),u(u_in){
  mode=1;
  max_dbl = 1e+20;//std::numeric_limits<double >::max();
  min_dbl = -1e+20;//std::numeric_limits<double >::min();
  num_eq=0; num_ineq=0;
  nnz_eq=0; nnz_ineq=0;
  problem_name = "noname";
  ql_wanted = 0;
  A_d = B_d = NULL;

 }

 ~QPFormatConverter(){
  if(mode==1 || mode==3){
   delete []a_eq;
   delete []a_ineq;
   allocateAC(H,0,0,0,FALSE);
   allocateAC(A_eq,0,0,0,FALSE);
   allocateAC(A_ineq,0,0,0,FALSE);
   allocateAC(AT_eq,0,0,0,FALSE);
   allocateAC(AT_ineq,0,0,0,FALSE);
   allocateAC(AT,0,0,0,FALSE);
  }
  if(mode==2 || mode==3){
   delete []l;
   delete []u;
   allocateAC(A,0,0,0,FALSE);
   allocateAC(H_full,0,0,0,FALSE);
  }
  if(mode == 0){
   allocateAC(H,0,0,0,FALSE);
   if(num_eq>0){
    allocateAC(A_eq,0,0,0,FALSE);
    delete []a_eq;
    allocateAC(AT_eq,0,0,0,FALSE);
    if(ql_wanted){
     allocateAC(H_full,0,0,0,FALSE);
    }
   } else {
    delete AT_eq;
    delete A_eq;
   }
   if(num_ineq>0){
    allocateAC(A_ineq,0,0,0,FALSE);
    delete []a_ineq;
    allocateAC(AT_ineq,0,0,0,FALSE);
   } else {
    delete A_ineq;
    delete AT_ineq;
   }
  }
  delete []q;
  if(ql_wanted){
   delete []AB_d;
   delete []ab_eqineq;
   delete []H_d;
   delete []A_d;
   delete []B_d;
  }
 }

 int read_IE_format(std::string hessian_file,
                    std::string linear_file,
                    std::string eq_file,
                    std::string eq_file_u,
                    std::string ineq_file,
                    std::string ineq_file_u){
  problem_name = linear_file;
  H = new CSC;
  CSC *H_tmp = new CSC;
  A_eq = new CSC;
  A_ineq = new CSC;
  AT_eq = new CSC;
  AT_ineq = new CSC;
  if (!readMatrix(hessian_file,H_tmp->ncol,H_tmp->nzmax,H_tmp->p,
                  H_tmp->i,H_tmp->x))
   return 0;
  //print_csc("HH\n ",H_tmp->ncol,H_tmp->p,H_tmp->i,H_tmp->x);
  bool is_expanded = expandMatrix(H_tmp->ncol,H_tmp->nzmax,H_tmp->p,
                                  H_tmp->i,H_tmp->x,
                                  H->nzmax,H->p,
                                  H->i,H->x);
  if(is_expanded){
   //H has the expanded version already.
   H->ncol = H->nrow = H_tmp->ncol;
   allocateAC(H_tmp,0,0,0,FALSE);
  } else{
   H->x = H_tmp->x;
   H->p = H_tmp->p;
   H->i = H_tmp->i;
   H->nzmax = H_tmp->nzmax;
   H->ncol = H->nrow = H_tmp->ncol;
  }
  H->nzmax = H->p[H->ncol];
  H->stype=-1;
  H->packed=1;
  //TODO: if it is full symmetric, make it half
/*  CSC *lower_H = computeLowerTriangular(H);
  if(lower_H){
   allocateAC(H_tmp,0,0,0,FALSE);
   H->ncol = H->nrow = lower_H->ncol;
   H->x = lower_H->x;
   H->p = lower_H->p;
   H->i = lower_H->i;
   H->nzmax = lower_H->nzmax;
  }*/

  if(eq_file != "none"){
   if (!readMatrix_rect(eq_file,A_eq->nrow,A_eq->ncol,A_eq->nzmax,
                        A_eq->p,A_eq->i,A_eq->x))
    return 0;
   if(A_eq->nrow>0){
    A_eq->nzmax = A_eq->p[A_eq->ncol];
    num_eq = A_eq->nrow;
    a_eq = new double[A_eq->nrow];
    read_vector(eq_file_u,A_eq->nrow,a_eq);
    transpose_unsym(A_eq->nrow,A_eq->ncol,A_eq->p,A_eq->i,
                    A_eq->x,AT_eq->nrow,AT_eq->ncol,AT_eq->p,AT_eq->i,
                    AT_eq->x);
   }else{
    A_eq->nzmax = A_eq->nrow = A_eq->ncol = 0;
    A_eq->p = A_eq->i = NULL;
    A_eq->x = a_eq = NULL;
    num_eq = 0;
   }
  }else{
   A_eq->nzmax = A_eq->nrow = A_eq->ncol = 0;
   A_eq->p = A_eq->i = NULL;
   A_eq->x = a_eq = NULL;
   num_eq = 0;
  }

  if(ineq_file != "none"){
   if (!readMatrix_rect(ineq_file,A_ineq->nrow,A_ineq->ncol,A_ineq->nzmax,
     A_ineq->p,A_ineq->i,A_ineq->x))
    return 0;
   if(A_ineq->nrow>0){
    A_ineq->nzmax = A_ineq->p[A_ineq->ncol];
    num_ineq = A_ineq->nrow;

    a_ineq = new double[A_ineq->nrow];
    read_vector(ineq_file_u,A_ineq->nrow,a_ineq);
    transpose_unsym(A_ineq->nrow,A_ineq->ncol,A_ineq->p,A_ineq->i,
                    A_ineq->x,AT_ineq->nrow,AT_ineq->ncol,AT_ineq->p,
                    AT_ineq->i,AT_ineq->x);
   }else{
    A_ineq->nzmax = A_ineq->nrow = A_ineq->ncol = 0;
    A_ineq->p = A_ineq->i = NULL;
    A_ineq->x = a_ineq = NULL;
    num_ineq = 0;
   }
  }else{
   A_ineq->nzmax = A_ineq->nrow = A_ineq->ncol = 0;
   A_ineq->p = A_ineq->i = NULL;
   A_ineq->x = a_ineq = NULL;
   num_ineq = 0;
  }

 q = new double[H->ncol];
 read_vector(linear_file,H->ncol,q);

  mode=0;
  delete H_tmp;
  return 1;
 }

 int read_bounded_format(std::string hessian_file,
                         std::string linear_file,
                         std::string ineq_file_l,
                         std::string ineq_file,
                         std::string ineq_file_u){
  problem_name = linear_file;
  H_full = new CSC;
  CSC *H_tmp = new CSC;
  A = new CSC;
  if (!readMatrix(hessian_file,H_tmp->ncol,H_tmp->nzmax,H_tmp->p,
                  H_tmp->i,H_tmp->x))
   return 0;
  //print_csc("HH\n ",H_tmp->ncol,H_tmp->p,H_tmp->i,H_tmp->x);
  bool is_expanded = expandMatrix(H_tmp->ncol,H_tmp->nzmax,H_tmp->p,
                                  H_tmp->i,H_tmp->x,
                                  H_full->nzmax,H_full->p,
                                  H_full->i,H_full->x);
  if(is_expanded){
   H_full->ncol = H_full->nrow = H_tmp->ncol;
   allocateAC(H_tmp,0,0,0,FALSE);
  } else{
   H_full->x = H_tmp->x;
   H_full->p = H_tmp->p;
   H_full->i = H_tmp->i;
   H_full->nzmax = H_tmp->nzmax;
   H_full->ncol = H_full->nrow = H_tmp->ncol;
   delete H_tmp;
  }

  if (!readMatrix_rect(ineq_file,A->nrow,A->ncol,A->nzmax,A->p,A->i,A->x))
   return 0;
  H_full->nzmax = H_full->p[H_full->ncol];
  A->nzmax = A->p[A->ncol];
  q = new double[H_full->ncol];
  read_vector(linear_file,H_full->ncol,q);
  l = new double[A->nrow];
  read_vector(ineq_file_l,A->nrow,l);
  u = new double[A->nrow];
  read_vector(ineq_file_u,A->nrow,u);
  mode=3;
  return 1;
 }

 int print_bounded_format(){
  print_csc("\nObjective function: \n",H_full->ncol,H_full->p,H_full->i,
    H_full->x);
  print_csc("\nConstraint matrix: \n",A->ncol,A->p,A->i,A->x);
  print_vec("\nLower bounds: \n",0,A->nrow,l);
  print_vec("\nUpper bounds: \n",0,A->nrow,u);
 }

 int B2IE(){
  if(mode!=1 && mode!=3)
   return 0;
  // Objective function
  H = new CSC;
  make_lower(H_full->ncol,H_full->nzmax,H_full->p,H_full->i,H_full->x,
             H->ncol,H->nzmax,H->p,H->i,H->x);

  // Constraints
  AT = new CSC;
  transpose_unsym(A->nrow,A->ncol,A->p,A->i,A->x,AT->nrow,AT->ncol,
    AT->p,AT->i,AT->x);
  //print_csc("\nAT: \n",AT->ncol,AT->p,AT->i,AT->x);
  A_eq = new CSC;
  A_ineq = new CSC;
  AT_eq = new CSC;
  AT_ineq = new CSC;
  std::vector<int> eq_idx;
  std::vector<constraint*> ineq_dx;
  int h_dim = H_full->ncol;
  int *col_cnt_A_eq = new int[h_dim]();
  int *col_cnt_A_ineq = new int[h_dim]();
  a_eq = new double[2*A->nrow];//FIXME: allocate it smarter
  a_ineq = new double[2*A->nrow];

  for (int i = 0; i < A->nrow; ++i) {
/*   l[i] = std::abs(l[i]) < 1e-14 ? 0 : l[i];
   u[i] = std::abs(u[i]) < 1e-14 ? 0 : u[i];*/
   // Invalid constraint
   if ((is_equal(l[i], min_dbl) && is_equal(u[i], max_dbl)) ||
       (is_equal(l[i], max_dbl) && is_equal(u[i], max_dbl)) ||
       (is_equal(l[i], min_dbl) && is_equal(u[i], min_dbl))) {
    continue;
   }
   if (is_equal(l[i], u[i])) {//eq

    a_eq[num_eq] = l[i] ;
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
  assert(num_ineq <= 2*A->nrow);
  assert(num_eq <= 2*A->nrow);
  // building the equality constraint matrix
  AT_eq->nrow=h_dim;
  AT_eq->ncol=num_eq;
  AT_eq->p = new int[num_eq+1];
  AT_eq->i = new int[nnz_eq];
  AT_eq->x = new double[nnz_eq];
  assert(eq_idx.size() == num_eq);
  AT_eq->p[0]=0;
  for (int ll = 0; ll < eq_idx.size(); ++ll) {
   int cur_idx = eq_idx[ll];
   int nnz_cur_row = AT->p[cur_idx+1] - AT->p[cur_idx];
   AT_eq->p[ll+1] = AT_eq->p[ll] + nnz_cur_row;
   for (int ii = AT_eq->p[ll], jj=AT->p[cur_idx];
        ii < AT_eq->p[ll + 1]; ++ii, ++jj) {
    AT_eq->i[ii] = AT->i[jj];
    AT_eq->x[ii] = AT->x[jj];
   }
  }
  transpose_unsym(AT_eq->nrow,AT_eq->ncol,AT_eq->p,AT_eq->i,
                  AT_eq->x,A_eq->nrow,A_eq->ncol,A_eq->p,A_eq->i,A_eq->x);
  A_eq->nzmax = nnz_eq;
  //building the ineq constraint matrix
  AT_ineq->nrow=h_dim;
  AT_ineq->ncol=num_ineq;
  AT_ineq->p = new int[num_ineq+1];
  AT_ineq->i = new int[nnz_ineq];
  AT_ineq->x = new double[nnz_ineq];
  assert(ineq_dx.size() == num_ineq);
  AT_ineq->p[0]=0;
  for (int ll = 0; ll < ineq_dx.size(); ++ll) {
   int cur_idx = ineq_dx[ll]->idx_no;
   double cur_coef = ineq_dx[ll]->coef;
   int nnz_cur_row = AT->p[cur_idx+1] - AT->p[cur_idx];
   AT_ineq->p[ll+1] = AT_ineq->p[ll] + nnz_cur_row;
   for (int ii = AT_ineq->p[ll], jj=AT->p[cur_idx];
   ii < AT_ineq->p[ll + 1]; ++ii, ++jj) {
    AT_ineq->i[ii] = AT->i[jj];
    AT_ineq->x[ii] = cur_coef*AT->x[jj];
   }
  }

  transpose_unsym(AT_ineq->nrow,AT_ineq->ncol,AT_ineq->p,AT_ineq->i,
                  AT_ineq->x,A_ineq->nrow,A_ineq->ncol,A_ineq->p,
                  A_ineq->i,A_ineq->x);
  A_ineq->nzmax = nnz_ineq;
  A_eq->stype=0;A_ineq->stype=0;
  H_full->stype=0;
  eq_idx.clear();
  for (int ii = 0; ii < ineq_dx.size(); ++ii) {
   delete ineq_dx[ii];
  }
  delete []col_cnt_A_ineq;
  delete []col_cnt_A_eq;
 }

 int print_IE_format(){

  print_csc("\nObjective function: \n",H->ncol,H->p,H->i,
            H->x);
  print_vec("\nq: \n",0,H->ncol,q);
  print_csc("\nEq matrix: \n",A_eq->ncol,A_eq->p,A_eq->i,A_eq->x);
  print_csc("\nIneq matrix: \n",A_ineq->ncol,A_ineq->p,A_ineq->i,A_ineq->x);
  print_vec("\nEq bounds: \n",0,num_eq,a_eq);
  print_vec("\nIneq bounds: \n",0,num_ineq,a_ineq);
 }

 int IE_export_to_dense(std::string mat_name){

  if(mode!=1 && mode!=3 && mode!=0)
   return 0;

  int total_cnst = A_ineq->nrow + A_eq->nrow;
  if(ql_wanted){
   if(mode != 0) {
    AB_d = new double[(total_cnst) * H_full->ncol]();
    ab_eqineq = new double[total_cnst];
   }else {
    H_full = new CSC;
    int status=0;
    CSC *HT = ptranspose(H, 2, NULL, NULL, 0, status);
    int nnzH = H->nzmax, nnz_full=0;
    make_full(H->ncol,nnzH,H->p,H->i,H->x,HT->p,HT->i,HT->x,nnz_full,
      H_full->p,H_full->i,H_full->x);
    H_full->nrow=H_full->ncol=H->ncol;H_full->stype=0;H_full->packed=1;
    AB_d = new double[(total_cnst) * H_full->ncol]();
    ab_eqineq = new double[total_cnst];
    allocateAC(HT,0,0,0,FALSE);
   }
  }
  H_d = new double [H_full->ncol*H_full->ncol]();
  sparse2dense(H_full,H_d);
  //write_dense(mat_name+"_P.dns",H_full->nrow,H_full->ncol,H_d);
  //write_vector(mat_name+"_q.dns",H_full->ncol,q);
  if(A_eq->nrow > 0){
   A_d = new double[A_eq->nrow*A_eq->ncol]();
   sparse2dense(A_eq,A_d);
   if(ql_wanted){
    for (int i = 0; i < A_eq->ncol; ++i) {
     for (int j = 0; j < A_eq->nrow; ++j) {
      AB_d[i*total_cnst + j] = -A_d[i*A_eq->nrow+j];
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
  if(A_ineq->nrow > 0){
   B_d = new double[A_ineq->nrow*A_ineq->ncol]();
   sparse2dense(A_ineq,B_d);
   if(ql_wanted){
    for (int i = 0; i < A_ineq->ncol; ++i) {
     for (int j = 0; j < A_ineq->nrow; ++j) {
      AB_d[i*total_cnst + j + A_eq->nrow] = -B_d[i*A_ineq->nrow+j];
     }
    }
    for (int k = 0; k < A_ineq->nrow; ++k) {
     ab_eqineq[k+A_eq->nrow] = a_ineq[k];
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
  if(!ql_wanted)
   delete []H_d;
 }

 void print_log(){
  std::setprecision(48);
  if(mode==1 || mode==3 || mode==0){
   std::cout<<problem_name<<","<<H->ncol<<","<<H->nzmax<<",";
   std::cout<<A_eq->nrow<<","<<A_eq->nzmax<<","<<A_ineq->nrow<<",";
   std::cout<<A_ineq->nzmax<<",";
  }
 }

};
#endif //PARS_QP_FORMAT_CONVERTER_H
