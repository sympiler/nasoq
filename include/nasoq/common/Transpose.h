//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_TRANSPOSE_H
#define CHOLOPENMP_TRANSPOSE_H

#include <cstdint>
#include <cstdio>

#include "nasoq/common/def.h"
#include "nasoq/common/SparseUtils.h"

namespace nasoq {
/* ========================================================================== */
/* === t_cholmod_transpose_sym ============================================== */
/* ========================================================================== */

/* Compute F = A' or A (p,p)', where A is symmetric and F is already allocated.
 * The complex case performs either the array transpose or complex conjugate
 * transpose.
 *
 * workspace:  Iwork (nrow) if Perm NULL, Iwork (2*nrow) if Perm non-NULL.
 */

 int transpose_sym_real
   (
     /* ---- input ---- */
     CSC *A, /* matrix to transpose */
     int *Perm,  /* size n, if present (can be NULL) */
     /* ---- output --- */
     CSC *F, /* F = A' or A(p,p)' */
     /* --------------- */
     //  cholmod_common *Common
     int *Wi,
     int *Pinv,
     int &status
   );

/* ========================================================================== */
/* === cholmod_transpose_sym ================================================ */
/* ========================================================================== */

/* Compute F = A' or A (p,p)', where A is symmetric and F is already allocated.
 * See cholmod_transpose for a simpler routine.
 *
 * workspace:  Iwork (nrow) if Perm NULL, Iwork (2*nrow) if Perm non-NULL.
 */

 int transpose_sym
   (
     /* ---- input ---- */
     CSC *A, /* matrix to transpose */
     int values,  /* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values */
     int *Perm,  /* size nrow, if present (can be NULL) */
     /* ---- output --- */
     CSC *F, /* F = A' or A(p,p)' */
     /* --------------- */
     //cholmod_common *Common
     int status
   );




/* ========================================================================== */
/* === cholmod_ptranspose =================================================== */
/* ========================================================================== */

/* Return A' or A(p,p)' if A is symmetric.  Return A', A(:,f)', or A(p,f)' if
 * A is unsymmetric.
 *
 * workspace:
 * Iwork (MAX (nrow,ncol)) if unsymmetric and fset is non-NULL
 * Iwork (nrow) if unsymmetric and fset is NULL
 * Iwork (2*nrow) if symmetric and Perm is non-NULL.
 * Iwork (nrow) if symmetric and Perm is NULL.
 *
 * A simple worst-case upper bound on the workspace is nrow+ncol.
 */

 CSC *ptranspose
   (
     /* ---- input ---- */
     CSC *A, /* matrix to transpose */
     int values,  /* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values */
     int *Perm,  /* if non-NULL, F = A(p,f) or A(p,p) */
     int *fset,  /* subset of 0:(A->ncol)-1 */
     size_t fsize, /* size of fset */
     /* --------------- */
     //cholmod_common *Common
     int status
   );

/* ========================================================================== */
/* === cholmod_transpose ==================================================== */
/* ========================================================================== */

/* Returns A'.  See also cholmod_ptranspose below. */

 CSC *transposeC
   (
     /* ---- input ---- */
     CSC *A, /* matrix to transpose */
     int values,  /* 2: complex conj. transpose, 1: array transpose,
			   0: do not transpose the numerical values
			   (returns its result as CHOLMOD_PATTERN) */
     /* --------------- */
     //cholmod_common *Common
     int &status
   );


//***********CSPARSE rountines
 double cumsum(int *p, int *c, int n);

 int transpose(int n, int m, int *Ap, int *Ai, double *Ax, int values,
               int *Cp, int *Ci, double *Cx);


}
#endif //CHOLOPENMP_TRANSPOSE_H
