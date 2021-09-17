//
// Created by kazem on 3/4/18.
//


// The code in this file is obtained from Suitesparse.
#ifndef PROJECT_SPMV_CSC_H
#define PROJECT_SPMV_CSC_H

#include <cstddef>

namespace nasoq{
#define IS_NAN(x)	(((x) != (x)) || (((x) < (x))))
#define IS_ZERO_R(x)	(((x) == 0.) && !IS_NAN(x))
#define IS_NONZERO_R(x)	(((x) != 0.) || IS_NAN(x))
#define TEMPLATE(name)			r_ ## name
#define ASSEMBLE(x,z,p,ax,az,q)		x [p] += ax [q]
#define ASSIGN(x,z,p,ax,az,q)			x [p]  = ax [q]
#define ASSIGN_CONJ(x,z,p,ax,az,q)		x [p]  = ax [q]
#define ASSIGN_REAL(x,p,ax,q)			x [p]  = ax [q]
#define XTYPE_OK(type)			((type) == CHOLMOD_REAL)
#define IS_NONZERO(ax,az,q)			IS_NONZERO_R (ax [q])
#define IS_ZERO(ax,az,q)			IS_ZERO_R (ax [q])
#define IS_ONE(ax,az,q)			(ax [q] == 1)
#define MULT(x,z,p, ax,az,q, bx,bz,r)		x [p]  = ax [q] * bx [r]
#define MULTADD(x,z,p, ax,az,q, bx,bz,r)	x [p] += ax [q] * bx [r]
#define MULTSUB(x,z,p, ax,az,q, bx,bz,r)	x [p] -= ax [q] * bx [r]
#define MULTADDCONJ(x,z,p, ax,az,q, bx,bz,r)	x [p] += ax [q] * bx [r]
#define MULTSUBCONJ(x,z,p, ax,az,q, bx,bz,r)	x [p] -= ax [q] * bx [r]
#define ADD(x,z,p, ax,az,q, bx,bz,r)		x [p]  = ax [q] + bx [r]
#define ADD_REAL(x,p, ax,q, bx,r)		x [p]  = ax [q] + bx [r]
#define CLEAR(x,z,p)				x [p]  = 0
#define CLEAR_IMAG(x,z,p)
#define DIV(x,z,p,ax,az,q)			x [p] /= ax [q]
#define LLDOT(x,p, ax,az,q)			x [p] -= ax [q] * ax [q]
#define PRINT(k,x,z,p)			PRK(k, ("%24.16e", x [p]))
#define ADVANCE(x,z,d) x += d
#define DIV_REAL(x,z,p, ax,az,q, bx,r)	x [p] = ax [q] / bx [r]
#define MULT_REAL(x,z,p, ax,az,q, bx,r)	x [p] = ax [q] * bx [r]
#define LDLDOT(x,p, ax,az,q, bx,r)		x [p] -=(ax[q] * ax[q])/ bx[r]


int spmv_csc (int n, size_t *Ap, int *Ai, double *Ax,
              double *x, double *y);

int spmv_csc_small (int m, int n, int *Ap, int *Ai, double *Ax,
              double *x, double *y);


/*
 * alpha * A * X + beta * Y
 */

void spmv_csc_sym_one_int (
  /* ---- input ---- */
  int ncol,
  int *Ap,
  int *Ai,
  double *Ax,
  int stype,
  double alpha [2],   /* scale factor for A */
  double beta [2],    /* scale factor for Y */
  int xcol,
  double *X,	/* dense matrix to multiply */
  /* ---- in/out --- */
  double *Y	/* resulting dense matrix */
);


/*
 * alpha * A * X + beta * Y
 */

void spmv_csc_sym_one (
  /* ---- input ---- */
  int ncol,
  size_t *Ap,
  int *Ai,
  double *Ax,
  int stype,
  double alpha [2],   /* scale factor for A */
  double beta [2],    /* scale factor for Y */
  int xcol,
  double *X,	/* dense matrix to multiply */
  /* ---- in/out --- */
  double *Y	/* resulting dense matrix */
);


/*
 * alpha * A * X + beta * Y
 */

void spmv_csc_sym (
    /* ---- input ---- */
    int ncol,
    size_t *Ap,
    int *Ai,
    double *Ax,
    int stype,
    double alpha [2],   /* scale factor for A */
    double beta [2],    /* scale factor for Y */
    int xcol,
    double *X,	/* dense matrix to multiply */
    /* ---- in/out --- */
    double *Y	/* resulting dense matrix */
  );
}
#endif //PROJECT_SPMV_CSC_H
