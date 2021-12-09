//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/Norm.h"

#include "nasoq/matrixVector/spmv_CSC.h"

#include <cmath>
#include <algorithm>


namespace nasoq {

 double norm_dense(int nrow, int ncol, double *X, int norm) {
  if (norm < 0 || norm > 2 || (norm == 2 && ncol > 1)) {
   std::cout << "NORM error";
   return -1;
  }
  int d = nrow;
  double xnorm = 0, s = 0;
  /* infinity-norm = max row sum, using stride-1 access of X */
  if (norm == 0) {
   /* infinity-norm = max row sum, using stride-d access of X */
   for (int i = 0; i < nrow; i++) {
    s = 0;
    for (int j = 0; j < ncol; j++) {
     s += std::abs(X[i + j * d]);
    }
    if ((!std::isfinite(s) || s > xnorm) && std::isfinite(xnorm)) {
     xnorm = s;
    }
   }
  } else if (norm == 1) {
   /* 1-norm = max column sum */
   for (int j = 0; j < ncol; j++) {
    s = 0;
    for (int i = 0; i < nrow; i++) {
     s += std::abs(X[i + j * d]);
    }
    if ((!std::isfinite(s) || s > xnorm) && std::isfinite(xnorm)) {
     xnorm = s;
    }
   }
  } else {
   /* 2-norm = sqrt (sum (X.^2)) */
   for (int i = 0; i < nrow; i++) {
    double x = X[i];
    xnorm += x * x;
   }
   xnorm = sqrt(xnorm);
  }
  return (xnorm);
 }

 double norm_sparse(int ncol, size_t *Ap, int *Ai, double *Ax, int stype, int norm) {
  double anorm, s;
  int p = 0, pend, i = 0;
  bool packed = true;
  double *W = new double[ncol]();
  /* check inputs */
  if (norm < 0 || norm > 1) {
   std::cout << "invalid norm\n";
   return -1;
  }
  anorm = 0;
  if (stype > 0) {
   /* A is symmetric with upper triangular part stored */
   /* infinity-norm = 1-norm = max row/col sum */
   for (int j = 0; j < ncol; j++) {
    p = Ap[j];
    pend = Ap[j + 1];
    for (; p < pend; p++) {
     i = Ai[p];
     s = std::abs(Ax[p]);
     if (i == j) {
      W[i] += s;
     } else if (i < j) {
      W[i] += s;
      W[j] += s;
     }
    }
   }
  } else if (stype < 0) {
   /* A is symmetric with lower triangular part stored */
   /* infinity-norm = 1-norm = max row/col sum */
   for (int j = 0; j < ncol; j++) {
    p = Ap[j];
    pend = Ap[j + 1];
    for (; p < pend; p++) {
     i = Ai[p];
     s = std::abs(Ax[p]);
     if (i == j) {
      W[i] += s;
     } else if (i > j) {
      W[i] += s;
      W[j] += s;
     }
    }
   }
  }
  /* compute the max row sum */
  if (stype || norm == 0) {
   for (i = 0; i < ncol; i++) {
    s = W[i];
    if ((std::isfinite(s) || s > anorm) && std::isfinite(anorm)) {
     anorm = s;
    }
   }
  }
  delete[]W;
  return (anorm);
 }

 double norm_sparse_int(int ncol, int *Ap, int *Ai, double *Ax, int stype, int norm) {
  double anorm, s, frob = 0;
  int p = 0, pend, i = 0;
  bool packed = true;
  double *W = new double[ncol]();
  /* check inputs */
  if (norm < 0 || norm > 3) {
   std::cout << "invalid norm\n";
   return -1;
  }
  anorm = 0;
  if (stype > 0) {
   /* A is symmetric with upper triangular part stored */
   /* infinity-norm = 1-norm = max row/col sum */
   for (int j = 0; j < ncol; j++) {
    p = Ap[j];
    pend = Ap[j + 1];
    for (; p < pend; p++) {
     i = Ai[p];
     s = std::abs(Ax[p]);
     if (i == j) {
      W[i] += s;
      frob += (s * s);
     } else if (i < j) {
      W[i] += s;
      W[j] += s;
      frob += (2 * s * s);
     }
    }
   }
  } else if (stype < 0) {
   /* A is symmetric with lower triangular part stored */
   /* infinity-norm = 1-norm = max row/col sum */
   for (int j = 0; j < ncol; j++) {
    p = Ap[j];
    pend = Ap[j + 1];
    for (; p < pend; p++) {
     i = Ai[p];
     s = std::abs(Ax[p]);
     if (i == j) {
      W[i] += s;
      frob += (s * s);
     } else if (i > j) {
      W[i] += s;
      W[j] += s;
      frob += (2 * s * s);
     }
    }
   }
  } else if (stype == 0) {
   for (int j = 0; j < ncol; j++) {
    for (int k = Ap[j]; k < Ap[j + 1]; ++k) {
     p = Ai[k];
     s = std::abs(Ax[p]);
     W[p] += s;
     frob += s * s;
    }
   }
  }
  /* compute the max row sum */
  if (norm == 0) {
   for (i = 0; i < ncol; i++) {
    s = W[i];
    if ((!std::isfinite(s) || s > anorm) && std::isfinite(anorm)) {
     anorm = s;
    }
   }
  }
  if (norm == 2) {
   if (std::isfinite(s)) {
    anorm = sqrt(frob);
   }
  }
  delete[]W;
  return (anorm);
 }

 double error_L2(int n, double *expected, double *x) {
  double x_error = 0.0;
  double exp_sum = 0.0;
  for (int i = 0; i < n; i++) {
   x_error += pow(expected[i] - x[i], 2);
  }
  x_error = sqrt(x_error);
  return x_error;
 }

 double error_relative_L2(int n, double *expected, double *x) {
  double x_error = 0.0;
  double exp_sum = 0.0;
  for (int i = 0; i < n; i++) {
   x_error += pow(expected[i] - x[i], 2);
   exp_sum += pow(expected[i], 2);
  }
  x_error = sqrt(x_error / exp_sum);
  return x_error;
 }

 double error_max(int n, double *expected, double *x) {
  double x_error = 0.0;
  double exp_max = 0.0;
  for (int i = 0; i < n; i++) {
   if (!std::isfinite(x[i]) || !std::isfinite(expected[i]) ||
       expected[i] == 0)
    return -1;
   x_error = std::abs(expected[i] - x[i]) / std::abs(expected[i]);
   exp_max = x_error > exp_max ? x_error : exp_max;
  }
  return exp_max;
 }

 double precision(int n, size_t *Ap, int *Ai, double *Ax, double *x, double *b) {
  double precision = 0;
  double alp[2] = {1.0, 0};
  double bet[2] = {1.0, 0};
  double *tmp = new double[n];
  for (int i = 0; i < n; ++i) {
   tmp[i] = -b[i];
  }
  spmv_csc_sym(n, Ap, Ai, Ax, -1, alp, bet, 1, x, tmp);
  double normAxMinusb = norm_dense(n, 1, tmp, 2);
  double normb = norm_dense(n, 1, b, 2);
  precision = normAxMinusb / normb;
  delete[]tmp;
  return precision;
 }

 double residual_norm(int n, size_t *Ap, int *Ai, double *Ax, double *x, double *b) {
  double precision = 0;
  double alp[2] = {1.0, 0};
  double bet[2] = {1.0, 0};
  double *tmp = new double[n];
  for (int i = 0; i < n; ++i) {
   tmp[i] = -b[i];
  }
  spmv_csc_sym(n, Ap, Ai, Ax, -1, alp, bet, 1, x, tmp);
  double normAxMinusb = norm_dense(n, 1, tmp, 0);
  return normAxMinusb;
 }

 double cs_norm1(int n, int *Ap, int *Ai, double *Ax) {
  int p, j;
  double norm = 0, s;
  if (!Ap || !Ax) return (-1); /* check inputs */
  for (j = 0; j < n; j++) {
   for (s = 0, p = Ap[j]; p < Ap[j + 1]; p++) s += fabs(Ax[p]);
   norm = std::max<double>(norm, s);
  }
  return (norm);
 }
}