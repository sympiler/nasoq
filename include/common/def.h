//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_DEF_H
#define CHOLOPENMP_DEF_H

#include <cassert>
#include <cstddef>
#ifndef _MSC_VER
#include "sys/param.h"

#else
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

namespace nasoq {
#define MKL 1
#define METIS 1



#define CS_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w, j) (w [j] < 0)
#define CS_MARK(w, j) { w [j] = CS_FLIP (w [j]) ; }
#define HEAD(k, j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
#define  CS_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define EMPTY -1
#define TRUE 1
#define FALSE 0
#define ASSERT(expression) (assert (expression))

#define CHOLMOD_NATURAL 0 /* use natural ordering */
#define CHOLMOD_GIVEN 1  /* use given permutation */
#define CHOLMOD_AMD 2  /* use minimum degree (AMD) */
#define CHOLMOD_METIS 3  /* use METIS' nested dissection */
#define CHOLMOD_NESDIS 4 /* use CHOLMOD's version of nested dissection:*/
/* node bisector applied recursively, followed
 * by constrained minimum degree (CSYMAMD or
 * CCOLAMD) */
#define CHOLMOD_COLAMD 5 /* use AMD for A, COLAMD for A*A' */

/* POSTORDERED is not a method, but a result of natural ordering followed by a
 * weighted postorder.  It is used for L->ordering, not method [ ].ordering. */
#define CHOLMOD_POSTORDERED 6 /* natural ordering, postordered. */
#define CHOLMOD_MAXMETHODS 9 /* maximum number of different methods that */
/* cholmod_analyze can try. Must be >= 9. */

/* find the sign: -1 if x < 0, 1 if x > 0, zero otherwise.
 * Not safe for NaN's */
#define SIGN(x) (((x) < 0) ? (-1) : (((x) > 0) ? 1 : 0))


/* xtype defines the kind of numerical values used: */
#define CHOLMOD_PATTERN 0 /* pattern only, no numerical values */
#define CHOLMOD_REAL 1  /* a real matrix */
#define CHOLMOD_COMPLEX 2 /* a complex matrix (ANSI C99 compatible) */
#define CHOLMOD_ZOMPLEX 3 /* a complex matrix (MATLAB compatible) */
#define Int_max INT_MAX
#define INVISIBLE -2
#define ROOT -1
/* ========================================================================== */
/* === Core/cholmod_sparse ================================================== */
/* ========================================================================== */

/* A sparse matrix stored in compressed-column form. */

 typedef struct sympiler_CSC {
  size_t nrow; /* the matrix is nrow-by-ncol */
  size_t ncol;
  size_t nzmax; /* maximum number of entries in the matrix */

  /* pointers to int or SuiteSparse_long: */
  int *p;  /* p [0..ncol], the column pointers */
  int *i;  /* i [0..nzmax-1], the row indices */

  /* for unpacked matrices only: */
  int *nz;  /* nz [0..ncol-1], the # of nonzeros in each col.  In
			 * packed form, the nonzero pattern of column j is in
	* A->i [A->p [j] ... A->p [j+1]-1].  In unpacked form, column j is in
	* A->i [A->p [j] ... A->p [j]+A->nz[j]-1] instead.  In both cases, the
	* numerical values (if present) are in the corresponding locations in
	* the array x (or z if A->xtype is CHOLMOD_ZOMPLEX). */

  /* pointers to double or float: */
  double *x;  /* size nzmax or 2*nzmax, if present */
  double *z;  /* size nzmax, if present */

  int stype;  /* Describes what parts of the matrix are considered:
			 *
	* 0:  matrix is "unsymmetric": use both upper and lower triangular parts
	*     (the matrix may actually be symmetric in pattern and value, but
	*     both parts are explicitly stored and used).  May be square or
	*     rectangular.
	* >0: matrix is square and symmetric, use upper triangular part.
	*     Entries in the lower triangular part are ignored.
	* <0: matrix is square and symmetric, use lower triangular part.
	*     Entries in the upper triangular part are ignored.
	*
	* Note that stype>0 and stype<0 are different for cholmod_sparse and
	* cholmod_triplet.  See the cholmod_triplet data structure for more
	* details.
	*/

  int itype;  /* CHOLMOD_INT:     p, i, and nz are int.
			 * CHOLMOD_INTLONG: p is SuiteSparse_long,
                         *                  i and nz are int.
			 * CHOLMOD_LONG:    p, i, and nz are SuiteSparse_long */

  int xtype;  /* pattern, real, complex, or zomplex */
  int dtype;  /* x and z are double or float */
  int sorted; /* TRUE if columns are sorted, FALSE otherwise */
  int packed; /* TRUE if packed (nz ignored), FALSE if unpacked
			 * (nz is required) */
  ~sympiler_CSC() {}

 } CSC;

/* ========================================================================== */
/* === Core/cholmod_factor ================================================== */
/* ========================================================================== */

/* A symbolic and numeric factorization, either simplicial or supernodal.
 * In all cases, the row indices in the columns of L are kept sorted. */

 typedef struct Sympiler_BCSC {
  /* ---------------------------------------------------------------------- */
  /* for both simplicial and supernodal factorizations */
  /* ---------------------------------------------------------------------- */

  size_t n;  /* L is n-by-n */

  size_t minor; /* If the factorization failed, L->minor is the column
			 * at which it failed (in the range 0 to n-1).  A value
			 * of n means the factorization was successful or
			 * the matrix has not yet been factorized. */

  /* ---------------------------------------------------------------------- */
  /* symbolic ordering and analysis */
  /* ---------------------------------------------------------------------- */

  int *Perm; /* size n, permutation used */
  int *ColCount; /* size n, column counts for simplicial L */

  int *IPerm;       /* size n, inverse permutation.  Only created by
                         * cholmod_solve2 if Bset is used. */

  /* ---------------------------------------------------------------------- */
  /* simplicial factorization */
  /* ---------------------------------------------------------------------- */

  size_t nzmax; /* size of i and x */

  size_t *p;  /* p [0..ncol], the column pointers */
  int *p_s;
  int *i;  /* i [0..nzmax-1], the row indices */
  double *x;  /* x [0..nzmax-1], the numerical values */
  double *x_s;
  void *z;
  void *nz;  /* nz [0..ncol-1], the # of nonzeros in each column.
			 * i [p [j] ... p [j]+nz[j]-1] contains the row indices,
			 * and the numerical values are in the same locatins
			 * in x. The value of i [p [k]] is always k. */

  void *next; /* size ncol+2. next [j] is the next column in i/x */
  void *prev; /* size ncol+2. prev [j] is the prior column in i/x.
			 * head of the list is ncol+1, and the tail is ncol. */

  /* ---------------------------------------------------------------------- */
  /* supernodal factorization */
  /* ---------------------------------------------------------------------- */

  /* Note that L->x is shared with the simplicial data structure.  L->x has
   * size L->nzmax for a simplicial factor, and size L->xsize for a supernodal
   * factor. */

  size_t nsuper; /* number of supernodes */
  size_t ssize; /* size of s, integer part of supernodes */
  size_t xsize, xsize_s; /* size of x, real part of supernodes */
  size_t maxcsize; /* size of largest update matrix */
  size_t maxesize; /* max # of rows in supernodes, excl. triangular part */

  int *super; /* size nsuper+1, first col in each supernode */ //SUP2Col
  int *col2Sup;
  size_t *pi;  /* size nsuper+1, pointers to integer patterns */ //row indices
  size_t *i_ptr;  /* size nsuper+1, pointers to integer patterns */ //row indices
  int *px;  /* size nsuper+1, pointers to real parts */ //
  int *s;  /* size ssize, integer part of supernodes */

  int *sParent; //Assembly Tree
  int *Parent;
  /* ---------------------------------------------------------------------- */
  /* factorization type */
  /* ---------------------------------------------------------------------- */

  int ordering; /* ordering method used */

  int is_ll;  /* TRUE if LL', FALSE if LDL' */
  int is_super; /* TRUE if supernodal, FALSE if simplicial */
  int is_monotonic; /* TRUE if columns of L appear in order 0..n-1.
			 * Only applicable to simplicial numeric types. */

  int itype; /* The integer arrays are Perm, ColCount, p, i, nz,
                 * next, prev, super, pi, px, and s.  If itype is
		 * CHOLMOD_INT, all of these are int arrays.
		 * CHOLMOD_INTLONG: p, pi, px are SuiteSparse_long, others int.
		 * CHOLMOD_LONG:    all integer arrays are SuiteSparse_long. */
  int xtype; /* pattern, real, complex, or zomplex */
  int dtype; /* x and z double or float */

  int useGPU; /* Indicates the symbolic factorization supports
		 * GPU acceleration */
  ~Sympiler_BCSC() {}

 } BCSC;

 void allocateLC(BCSC *L, int sw);

 void allocateAC(CSC *A, int nrow, int nnz, int sytpe, int sw);

/* define pour l'affichage */
#define SCOTCH_STRAT_DIRECT                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                       \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                        "m{rat=0.8,"                                \
  ""                          "vert=100,"                               \
  ""                          "low=h{pass=10},"                         \
  ""                          "asc=f{bal=0.2}};,"                       \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"

}
#endif //CHOLOPENMP_DEF_H
