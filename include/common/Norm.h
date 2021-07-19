//
// Created by Kazem on 11/24/18.
//

#ifndef PROJECT_NORM_H
#define PROJECT_NORM_H

#include <cmath>
#include <iostream>
#include <algorithm>

namespace nasoq {
 double norm_dense
   (
     /* ---- input ---- */
     int nrow,
     int ncol,
     double *X, /* matrix to compute the norm of */
     int norm  /* type of norm: 0: inf. norm, 1: 1-norm, 2: 2-norm */
   );

/*
 * Computes norm 0 of a sparse matrix
 */
 double norm_sparse
   (
     /* ---- input ---- */
     int ncol,
     size_t *Ap,
     int *Ai,
     double *Ax,
     int stype,  // 1 is upper, -1 is lower triangular
     int norm /* type of norm: 0: inf. norm, 1: 1-norm */

   );

/*
 * Computes norm 0 of a sparse matrix
 */
 double norm_sparse_int
   (
     /* ---- input ---- */
     int ncol,
     int *Ap,
     int *Ai,
     double *Ax,
     int stype,  // 1 is upper, -1 is lower triangular
     int norm /* type of norm: 0: inf. norm, 1: 1-norm ,
 * 2: frob norm*/

   );


/*
 * Compute L2  error norm
 */
 double error_L2(int n, double *expected, double *x);


/*
 * Compute L2 relative error norm
 */
 double error_relative_L2(int n, double *expected, double *x);

/*
 * computes maximum error
 */
 double error_max(int n, double *expected, double *x);

/*
 * Computing precision ||Ax-b|| / ||b||
 */
 double precision(int n,
                  size_t *Ap,
                  int *Ai,
                  double *Ax,
                  double *x,
                  double *b);

/*
 * inifinity norm of residual or ||Ax-b||inf
 */
 double residual_norm(int n,
                      size_t *Ap,
                      int *Ai,
                      double *Ax,
                      double *x,
                      double *b);


 double cs_norm1(int n, int *Ap, int *Ai, double *Ax);
}

#endif //PROJECT_NORM_H
