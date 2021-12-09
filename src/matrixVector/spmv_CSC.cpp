//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/matrixVector/spmv_CSC.h"

#include <algorithm>
#include <cassert>

namespace nasoq {

 int spmv_csc(int n, size_t *Ap, int *Ai, double *Ax, double *x, double *y) {
  int p, j;
  std::fill_n(y,n,0);
  if (!Ap || !x || !y) return (0) ;       /* check inputs */
  for (j = 0 ; j < n ; j++)
  {
   for (p = Ap [j] ; p < Ap [j+1] ; p++)
   {
    y [Ai [p]] += Ax [p] * x [j] ;
   }
  }
  return (1) ;
 }

 void spmv_csc_sym_one_int(int ncol, int *Ap, int *Ai, double *Ax, int stype, double *alpha, double *beta, int xcol,
                           double *X, double *Y) {

  double yx [8], xx [8], ax [2] ;
  double *Xz, *Yz,  *Wz ; //Xz and Yz are useless here
  size_t nx, ny, dx, dy ;
  int j, k, p, pend, kcol, i ;


  ny = nx = ncol;
  //double *w = new double[4*nx]();
  int nrow = ncol ;// only support square matrices
  kcol = xcol ;
  dy = xcol ;//FIXME: assuming stride is the same size
  dx = xcol ;

  /* Y = beta * Y */
  if (IS_ZERO (beta, betaz, 0)) {
   for (k = 0 ; k < kcol ; k++) {
    for (i = 0 ; i < ((int) ny) ; i++) {
     /* y [i] = 0. ; */
     CLEAR (Y, Yz, i) ;
    }
   }
  }
  else if (!IS_ONE (beta, betaz, 0)) {
   for (k = 0 ; k < kcol ; k++) {
    for (i = 0 ; i < ((int) ny) ; i++) {
     /* y [i] *= beta [0] ; */
     MULT (Y,Yz,i, Y,Yz,i, beta,betaz, 0) ;
    }
   }
  }

  if (IS_ZERO (alpha, alphaz, 0)) {
   /* nothing else to do */
   return ;
  }
  k = 0 ;
  /* Y += alpha * (A or A') * x, symmetric case (upper/lower) */
  /* Only the upper/lower triangular part and the diagonal of A is used.
   * Since both x and y are written to in the innermost loop, this
   * code can experience cache bank conflicts if x is used directly.
   * Thus, a copy is made of x, four columns at a time, if x has
   * four or more columns.
   */
  for (j = 0 ; j < ncol ; j++) {
   /* yj = 0. ; */
   CLEAR (yx,yz,0) ;
   /* xj = alpha [0] * x [j] ; */
   MULT (xx,xz,0, alpha,alphaz,0, X,Xz,j) ;
   p = Ap [j] ;
   pend = Ap [j+1];
   for ( ; p < pend ; p++) {
    i = Ai [p] ;
    if (i == j) {
     /* y [i] += Ax [p] * xj ; */
     MULTADD (Y,Yz,i, Ax,Az,p, xx,xz,0) ;
    }
    else if ((stype > 0 && i < j) || (stype < 0 && i > j)) {
     /* aij = Ax [p] ; */
     ASSIGN (ax,az,0, Ax,Az,p) ;
     /* y [i] += aij * xj ; */
     /* yj    += aij * x [i] ; */
     MULTADD     (Y,Yz,i, ax,az,0, xx,xz,0) ;
     MULTADDCONJ (yx,yz,0, ax,az,0, X,Xz,i) ;
    }
   }
   /* y [j] += alpha [0] * yj ; */
   MULTADD (Y,Yz,j, alpha,alphaz,0, yx,yz,0) ;
  }
  /* y += dy ; */
  /* x += dx ; */
  // ADVANCE (Y,Yz,dy) ;
  //ADVANCE (X,Xz,dx) ;
  k++ ;
 }

 void spmv_csc_sym_one(int ncol, size_t *Ap, int *Ai, double *Ax, int stype, double *alpha, double *beta, int xcol,
                       double *X, double *Y) {

  double yx [8], xx [8], ax [2] ;
  double *Xz, *Yz,  *Wz ; //Xz and Yz are useless here
  size_t nx, ny, dx, dy ;
  int j, k, p, pend, kcol, i ;


  ny = nx = ncol;
  double *w = new double[4*nx]();
  int nrow = ncol ;// only support square matrices
  kcol = xcol ;
  dy = xcol ;//FIXME: assuming stride is the same size
  dx = xcol ;

  /* Y = beta * Y */
  if (IS_ZERO (beta, betaz, 0)) {
   for (k = 0 ; k < kcol ; k++) {
    for (i = 0 ; i < ((int) ny) ; i++) {
     /* y [i] = 0. ; */
     CLEAR (Y, Yz, i) ;
    }
   }
  }
  else if (!IS_ONE (beta, betaz, 0)) {
   for (k = 0 ; k < kcol ; k++) {
    for (i = 0 ; i < ((int) ny) ; i++) {
     /* y [i] *= beta [0] ; */
     MULT (Y,Yz,i, Y,Yz,i, beta,betaz, 0) ;
    }
   }
  }

  if (IS_ZERO (alpha, alphaz, 0)) {
   /* nothing else to do */
   return ;
  }
  k = 0 ;
  /* Y += alpha * (A or A') * x, symmetric case (upper/lower) */
  /* Only the upper/lower triangular part and the diagonal of A is used.
   * Since both x and y are written to in the innermost loop, this
   * code can experience cache bank conflicts if x is used directly.
   * Thus, a copy is made of x, four columns at a time, if x has
   * four or more columns.
   */
  for (j = 0 ; j < ncol ; j++) {
   /* yj = 0. ; */
   CLEAR (yx,yz,0) ;
   /* xj = alpha [0] * x [j] ; */
   MULT (xx,xz,0, alpha,alphaz,0, X,Xz,j) ;
   p = Ap [j] ;
   pend = Ap [j+1];
   for ( ; p < pend ; p++) {
    i = Ai [p] ;
    if (i == j) {
     /* y [i] += Ax [p] * xj ; */
     MULTADD (Y,Yz,i, Ax,Az,p, xx,xz,0) ;
    }
    else if ((stype > 0 && i < j) || (stype < 0 && i > j)) {
     /* aij = Ax [p] ; */
     ASSIGN (ax,az,0, Ax,Az,p) ;
     /* y [i] += aij * xj ; */
     /* yj    += aij * x [i] ; */
     MULTADD     (Y,Yz,i, ax,az,0, xx,xz,0) ;
     MULTADDCONJ (yx,yz,0, ax,az,0, X,Xz,i) ;
    }
   }
   /* y [j] += alpha [0] * yj ; */
   MULTADD (Y,Yz,j, alpha,alphaz,0, yx,yz,0) ;
  }
  /* y += dy ; */
  /* x += dx ; */
  // ADVANCE (Y,Yz,dy) ;
  //ADVANCE (X,Xz,dx) ;
  k++ ;
 }

 void
 spmv_csc_sym(int ncol, size_t *Ap, int *Ai, double *Ax, int stype, double *alpha, double *beta, int xcol, double *X,
              double *Y) {

  double yx [8], xx [8], ax [2] ;
  double *Xz, *Yz,  *Wz ; //Xz and Yz are useless here
  size_t nx, ny, dx, dy ;
  int j, k, p, pend, kcol, i ;


  ny = nx = ncol;
  double *w = new double[4*nx]();
  int nrow = ncol ;// only support square matrices
  kcol = xcol ;
  dy = xcol ;//FIXME: assuming stride is the same size
  dx = xcol ;

  /* Y = beta * Y */
  if (IS_ZERO (beta, betaz, 0)) {
   for (k = 0 ; k < kcol ; k++) {
    for (i = 0 ; i < ((int) ny) ; i++) {
     /* y [i] = 0. ; */
     CLEAR (Y, Yz, i) ;
    }
    /* y += dy ; */
    //  ADVANCE (Y,Yz,dy) ;
   }
  }
  else if (!IS_ONE (beta, betaz, 0)) {
   for (k = 0 ; k < kcol ; k++) {
    for (i = 0 ; i < ((int) ny) ; i++) {
     /* y [i] *= beta [0] ; */
     MULT (Y,Yz,i, Y,Yz,i, beta,betaz, 0) ;
    }
    /* y += dy ; */
    ADVANCE (Y,Yz,dy) ;
   }
  }

  if (IS_ZERO (alpha, alphaz, 0)) {
   /* nothing else to do */
   return ;
  }
  k = 0 ;
  /* Y += alpha * (A or A') * x, symmetric case (upper/lower) */
  /* Only the upper/lower triangular part and the diagonal of A is used.
   * Since both x and y are written to in the innermost loop, this
   * code can experience cache bank conflicts if x is used directly.
   * Thus, a copy is made of x, four columns at a time, if x has
   * four or more columns.
   */
  if (kcol % 4 == 1) {
   for (j = 0 ; j < ncol ; j++) {
    /* yj = 0. ; */
    CLEAR (yx,yz,0) ;
    /* xj = alpha [0] * x [j] ; */
    MULT (xx,xz,0, alpha,alphaz,0, X,Xz,j) ;
    p = Ap [j] ;
    pend = Ap [j+1];
    for ( ; p < pend ; p++) {
     i = Ai [p] ;
     if (i == j) {
      /* y [i] += Ax [p] * xj ; */
      MULTADD (Y,Yz,i, Ax,Az,p, xx,xz,0) ;
     }
     else if ((stype > 0 && i < j) || (stype < 0 && i > j)) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i] += aij * xj ; */
      /* yj    += aij * x [i] ; */
      MULTADD     (Y,Yz,i, ax,az,0, xx,xz,0) ;
      MULTADDCONJ (yx,yz,0, ax,az,0, X,Xz,i) ;
     }
    }
    /* y [j] += alpha [0] * yj ; */
    MULTADD (Y,Yz,j, alpha,alphaz,0, yx,yz,0) ;
   }
   /* y += dy ; */
   /* x += dx ; */
   // ADVANCE (Y,Yz,dy) ;
   ADVANCE (X,Xz,dx) ;
   k++ ;
  }
  else if (kcol % 4 == 2) {
   for (j = 0 ; j < ncol ; j++) {
    /* yj0 = 0. ; */
    /* yj1 = 0. ; */
    CLEAR (yx,yz,0) ;
    CLEAR (yx,yz,1) ;
    /* xj0 = alpha [0] * x [j   ] ; */
    /* xj1 = alpha [0] * x [j+dx] ; */
    MULT (xx,xz,0, alpha,alphaz,0, X,Xz,j) ;
    MULT (xx,xz,1, alpha,alphaz,0, X,Xz,j+dx) ;
    p = Ap [j] ;
    pend = Ap [j+1] ;
    for ( ; p < pend ; p++) {
     i = Ai [p] ;
     if (i == j) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i   ] += aij * xj0 ; */
      /* y [i+dy] += aij * xj1 ; */
      MULTADD (Y,Yz,i,    ax,az,0, xx,xz,0) ;
      MULTADD (Y,Yz,i+dy, ax,az,0, xx,xz,1) ;
     }
     else if ((stype > 0 && i < j) || (stype < 0 && i > j)) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i   ] += aij * xj0 ; */
      /* y [i+dy] += aij * xj1 ; */
      /* yj0 += aij * x [i   ] ; */
      /* yj1 += aij * x [i+dx] ; */
      MULTADD     (Y,Yz,i,    ax,az,0, xx,xz,0) ;
      MULTADD     (Y,Yz,i+dy, ax,az,0, xx,xz,1) ;
      MULTADDCONJ (yx,yz,0,    ax,az,0, X,Xz,i) ;
      MULTADDCONJ (yx,yz,1,    ax,az,0, X,Xz,i+dx) ;
     }
    }
    /* y [j   ] += alpha [0] * yj0 ; */
    /* y [j+dy] += alpha [0] * yj1 ; */
    MULTADD (Y,Yz,j,    alpha,alphaz,0, yx,yz,0) ;
    MULTADD (Y,Yz,j+dy, alpha,alphaz,0, yx,yz,1) ;
   }
   /* y += 2*dy ; */
   /* x += 2*dx ; */
   ADVANCE (Y,Yz,2*dy) ;
   ADVANCE (X,Xz,2*dx) ;
   k += 2 ;
  }
  else if (kcol % 4 == 3) {
   for (j = 0 ; j < ncol ; j++) {
    /* yj0 = 0. ; */
    /* yj1 = 0. ; */
    /* yj2 = 0. ; */
    CLEAR (yx,yz,0) ;
    CLEAR (yx,yz,1) ;
    CLEAR (yx,yz,2) ;
    /* xj0 = alpha [0] * x [j     ] ; */
    /* xj1 = alpha [0] * x [j+  dx] ; */
    /* xj2 = alpha [0] * x [j+2*dx] ; */
    MULT (xx,xz,0, alpha,alphaz,0, X,Xz,j) ;
    MULT (xx,xz,1, alpha,alphaz,0, X,Xz,j+dx) ;
    MULT (xx,xz,2, alpha,alphaz,0, X,Xz,j+2*dx) ;
    p = Ap [j] ;
    pend = Ap [j+1] ;
    for ( ; p < pend ; p++) {
     i = Ai [p] ;
     if (i == j) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i     ] += aij * xj0 ; */
      /* y [i+  dy] += aij * xj1 ; */
      /* y [i+2*dy] += aij * xj2 ; */
      MULTADD (Y,Yz,i,      ax,az,0, xx,xz,0) ;
      MULTADD (Y,Yz,i+dy,   ax,az,0, xx,xz,1) ;
      MULTADD (Y,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
     }
     else if ((stype > 0 && i < j) || (stype < 0 && i > j)) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i     ] += aij * xj0 ; */
      /* y [i+  dy] += aij * xj1 ; */
      /* y [i+2*dy] += aij * xj2 ; */
      /* yj0 += aij * x [i     ] ; */
      /* yj1 += aij * x [i+  dx] ; */
      /* yj2 += aij * x [i+2*dx] ; */
      MULTADD     (Y,Yz,i,      ax,az,0, xx,xz,0) ;
      MULTADD     (Y,Yz,i+dy,   ax,az,0, xx,xz,1) ;
      MULTADD     (Y,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
      MULTADDCONJ (yx,yz,0,      ax,az,0, X,Xz,i) ;
      MULTADDCONJ (yx,yz,1,      ax,az,0, X,Xz,i+dx) ;
      MULTADDCONJ (yx,yz,2,      ax,az,0, X,Xz,i+2*dx) ;
     }
    }
    /* y [j     ] += alpha [0] * yj0 ; */
    /* y [j+  dy] += alpha [0] * yj1 ; */
    /* y [j+2*dy] += alpha [0] * yj2 ; */
    MULTADD (Y,Yz,j,      alpha,alphaz,0, yx,yz,0) ;
    MULTADD (Y,Yz,j+dy,   alpha,alphaz,0, yx,yz,1) ;
    MULTADD (Y,Yz,j+2*dy, alpha,alphaz,0, yx,yz,2) ;
   }
   /* y += 3*dy ; */
   /* x += 3*dx ; */
   ADVANCE (Y,Yz,3*dy) ;
   ADVANCE (X,Xz,3*dx) ;
   k += 3 ;
  }
  /* copy four columns of X into W, and put in row form */
  for ( ; k < kcol ; k += 4) {
   for (j = 0 ; j < ncol ; j++) {
    /* w [4*j  ] = x [j     ] ; */
    /* w [4*j+1] = x [j+  dx] ; */
    /* w [4*j+2] = x [j+2*dx] ; */
    /* w [4*j+3] = x [j+3*dx] ; */
    ASSIGN (w,Wz,4*j  , X,Xz,j     ) ;
    ASSIGN (w,Wz,4*j+1, X,Xz,j+dx  ) ;
    ASSIGN (w,Wz,4*j+2, X,Xz,j+2*dx) ;
    ASSIGN (w,Wz,4*j+3, X,Xz,j+3*dx) ;
   }
   for (j = 0 ; j < ncol ; j++) {
    /* yj0 = 0. ; */
    /* yj1 = 0. ; */
    /* yj2 = 0. ; */
    /* yj3 = 0. ; */
    CLEAR (yx,yz,0) ;
    CLEAR (yx,yz,1) ;
    CLEAR (yx,yz,2) ;
    CLEAR (yx,yz,3) ;
    /* xj0 = alpha [0] * w [4*j  ] ; */
    /* xj1 = alpha [0] * w [4*j+1] ; */
    /* xj2 = alpha [0] * w [4*j+2] ; */
    /* xj3 = alpha [0] * w [4*j+3] ; */
    MULT (xx,xz,0, alpha,alphaz,0, w,Wz,4*j) ;
    MULT (xx,xz,1, alpha,alphaz,0, w,Wz,4*j+1) ;
    MULT (xx,xz,2, alpha,alphaz,0, w,Wz,4*j+2) ;
    MULT (xx,xz,3, alpha,alphaz,0, w,Wz,4*j+3) ;
    p = Ap [j] ;
    pend = Ap [j+1] ;
    for ( ; p < pend ; p++) {
     i = Ai [p] ;
     assert(i < ncol);
     if (i == j) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i     ] += aij * xj0 ; */
      /* y [i+  dy] += aij * xj1 ; */
      /* y [i+2*dy] += aij * xj2 ; */
      /* y [i+3*dy] += aij * xj3 ; */
      MULTADD (Y,Yz,i     , ax,az,0, xx,xz,0) ;
      MULTADD (Y,Yz,i+dy  , ax,az,0, xx,xz,1) ;
      MULTADD (Y,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
      MULTADD (Y,Yz,i+3*dy, ax,az,0, xx,xz,3) ;
     }
     else if ((stype > 0 && i < j) || (stype < 0 && i > j)) {
      /* aij = Ax [p] ; */
      ASSIGN (ax,az,0, Ax,Az,p) ;
      /* y [i     ] += aij * xj0 ; */
      /* y [i+  dy] += aij * xj1 ; */
      /* y [i+2*dy] += aij * xj2 ; */
      /* y [i+3*dy] += aij * xj3 ; */
      /* yj0 += aij * w [4*i  ] ; */
      /* yj1 += aij * w [4*i+1] ; */
      /* yj2 += aij * w [4*i+2] ; */
      /* yj3 += aij * w [4*i+3] ; */
      MULTADD     (Y,Yz,i,      ax,az,0, xx,xz,0) ;
      MULTADD     (Y,Yz,i+dy,   ax,az,0, xx,xz,1) ;
      MULTADD     (Y,Yz,i+2*dy, ax,az,0, xx,xz,2) ;
      MULTADD     (Y,Yz,i+3*dy, ax,az,0, xx,xz,3) ;
      MULTADDCONJ (yx,yz,0,     ax,az,0, w,Wz,4*i) ;
      MULTADDCONJ (yx,yz,1,     ax,az,0, w,Wz,4*i+1) ;
      MULTADDCONJ (yx,yz,2,     ax,az,0, w,Wz,4*i+2) ;
      MULTADDCONJ (yx,yz,3,     ax,az,0, w,Wz,4*i+3) ;
     }
    }
    /* y [j     ] += alpha [0] * yj0 ; */
    /* y [j+  dy] += alpha [0] * yj1 ; */
    /* y [j+2*dy] += alpha [0] * yj2 ; */
    /* y [j+3*dy] += alpha [0] * yj3 ; */
    MULTADD (Y,Yz,j     , alpha,alphaz,0, yx,yz,0) ;
    MULTADD (Y,Yz,j+dy  , alpha,alphaz,0, yx,yz,1) ;
    MULTADD (Y,Yz,j+2*dy, alpha,alphaz,0, yx,yz,2) ;
    MULTADD (Y,Yz,j+3*dy, alpha,alphaz,0, yx,yz,3) ;
   }
   /* y += 4*dy ; */
   /* x += 4*dx ; */
   ADVANCE (Y,Yz,4*dy) ;
   ADVANCE (X,Xz,4*dx) ;
  }
  delete []w;
 }

 int spmv_csc_small(int m, int n, int *Ap, int *Ai, double *Ax, double *x, double *y) {
  int p, j;
  std::fill_n(y,m,0);
  if (!Ap || !x || !y) return (0) ;       /* check inputs */
  for (j = 0 ; j < n ; j++) {
   for (p = Ap [j] ; p < Ap [j+1] ; p++) {
    y [Ai [p]] += Ax [p] * x [j] ;
   }
  }
  return (1) ;
 }
}