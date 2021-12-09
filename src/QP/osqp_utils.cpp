//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/QP/osqp_utils.h"

#include <algorithm>
#include <cmath>

namespace nasoq {

 void vec_mult_scalar(double *a, double sc, int n) {
  int i;

  for (i = 0; i < n; i++) {
   a[i] *= sc;
  }
 }

 double vec_mean(const double *a, int n) {
  double mean = 0.0;
  int i;

  for (i = 0; i < n; i++) {
   mean += a[i];
  }
  mean /= (double) n;

  return mean;
 }

 void mat_inf_norm_cols_sym_triu(const CSC *M, double *E) {
  int i, j, ptr;
  double abs_x;

  // Initialize zero max elements
  for (j = 0; j < M->ncol; j++) {
   E[j] = 0.;
  }

  // Compute maximum across columns
  // Note that element (i, j) contributes to
  // -> Column j (as expected in any matrices)
  // -> Column i (which is equal to row i for symmetric matrices)
  for (j = 0; j < M->ncol; j++) {
   for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
    i = M->i[ptr];
    abs_x = std::abs(M->x[ptr]);
    E[j] = std::max(abs_x, E[j]);

    if (i != j) {
     E[i] = std::max(abs_x, E[i]);
    }
   }
  }
 }

 void mat_premult_diag(CSC *A, const double *d) {
  int j, i;

  for (j = 0; j < A->ncol; j++) {                // Cycle over columns
   for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
    A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
    // of d for row i
   }
  }
 }

 void mat_postmult_diag(CSC *A, const double *d) {
  int j, i;

  for (j = 0; j < A->ncol; j++) {                // Cycle over columns j
   for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
    A->x[i] *= d[j];                        // Scale by corresponding element
    // of d for column j
   }
  }
 }

 void vec_ew_prod(const double *a, const double *b, double *c, int n) {
  int i;

  for (i = 0; i < n; i++) {
   c[i] = b[i] * a[i];
  }
 }

 void vec_ew_recipr(const double *a, double *b, int n) {
  int i;

  for (i = 0; i < n; i++) {
   b[i] = (double) 1.0 / a[i];
  }
 }

 void vec_ew_sqrt(double *a, int n) {
  int i;

  for (i = 0; i < n; i++) {
   a[i] = std::sqrt(a[i]);
  }
 }

 void vec_set_scalar(double *a, double sc, int n) {
  int i;

  for (i = 0; i < n; i++) {
   a[i] = sc;
  }
 }

 void limit_scaling(double *D, int n) {
  int i;
  for (i = 0; i < n; i++) {
   D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
   D[i] = D[i] > MAX_SCALING ? MAX_SCALING : D[i];
  }
 }

 void mat_inf_norm_cols(const CSC *M, double *E) {
  int j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->ncol; j++) {
   E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->ncol; j++) {
   for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
    E[j] = std::max(std::abs(M->x[ptr]), E[j]);
   }
  }
 }

 void vec_ew_max_vec(const double *a, const double *b, double *c, int n) {
  int i;

  for (i = 0; i < n; i++) {
   c[i] = std::max(a[i], b[i]);
  }
 }

 void mat_inf_norm_rows(const CSC *M, double *E) {
  int i, j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->nrow; j++) {
   E[j] = 0.;
  }

  // Compute maximum across rows
  for (j = 0; j < M->ncol; j++) {
   for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
    i = M->i[ptr];
    E[i] = std::max(std::abs(M->x[ptr]), E[i]);
   }
  }
 }

 double vec_norm_inf(const double *v, int l) {
  int i;
  double abs_v_i;
  double max = 0.0;

  for (i = 0; i < l; i++) {
   abs_v_i = std::abs(v[i]);

   if (abs_v_i > max) max = abs_v_i;
  }
  return max;
 }

 void mat_mult_scalar(CSC *A, double sc) {
  int i, nnzA;

  nnzA = A->p[A->ncol];

  for (i = 0; i < nnzA; i++) {
   A->x[i] *= sc;
  }
 }

 void
 compute_inf_norm_cols_KKT(const CSC *P, const CSC *A, CSC *B, double *D, double *D_temp_A, double *E, double *F, int n) {
  // First half
  //  [ P ]
  //  [ A ]
  //  [ B ]
  mat_inf_norm_cols_sym_triu(P, D);
  mat_inf_norm_cols(A, D_temp_A);
  vec_ew_max_vec(D, D_temp_A, D, n);
  mat_inf_norm_cols(B, D_temp_A);
  vec_ew_max_vec(D, D_temp_A, D, n);
  // Second half
  //  [ A']
  //  [ B']
  //  [ 0 ]
  mat_inf_norm_rows(A, E);
  mat_inf_norm_rows(B, F);
 }
}