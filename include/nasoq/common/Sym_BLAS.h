//
// Created by kazem on 11/29/18.
//

#ifndef PROJECT_BLASKERNELS_H
#define PROJECT_BLASKERNELS_H

#include <cstdint>
#include <vector>
#include <cstdlib>
#include "nasoq/common/Reach.h"

#ifdef OPENBLAS
/*    #ifdef OB_INTERNAL
    #include "lapacke.h"
    #include "cblas.h"
    //#include "common_interface.h"
    #else*/
    //#include "openblas/f77blas.h"
#ifdef NASOQ_USE_CLAPACK
    #include "nasoq/clapacke/clapacke.h"
#else
    #include "openblas/lapacke.h"
#endif
    #include "openblas/cblas.h"
   // #endif
#else
#include "mkl.h"
#include <mkl_blas.h>
#include <mkl_lapacke.h>
#endif
namespace nasoq {
#  define VEC_SCAL(n, a, x, u){               \
    int i; double *pt,*p=(x);                    \
    for((pt)=(p+n);(p)<(pt);)                   \
      *((p)++)*= (a);                           \
  }

#ifdef OPENBLAS
#define SYM_DGEMM dgemm_
#define SYM_DTRSM dtrsm_
#define SYM_DGEMV dgemv_
#define SYM_DSCAL dscal_
#define SET_BLAS_THREAD(t) (openblas_set_num_threads(t))
#else
#define SYM_DGEMM dgemm
#define SYM_DTRSM dtrsm
#define SYM_DGEMV dgemv
#define SYM_DSCAL dscal

#define SET_BLAS_THREAD(t) (MKL_Domain_Set_Num_Threads(t, MKL_DOMAIN_BLAS))
#endif


 void
 sym_sytrf(double *A, int n, const int stride, int *nbpivot,
           double critere);

 double dot(int n, double *a, double *b);


/*
 * swaps a and b
 */
 void swap_vector(int n, double *a, double *b, int lda);

 void swap_int(int &a, int &b);

/*
 * Performs reverse ordering after sytrf
 * assume cur has full lower triangular.
 */
 int reorder_after_sytrf(int n, double *a, int lda, int *ipiv,
                         int *perm, double *D, int lda_D,
                         int *swap_vec, int *ws);

/*
 * Shoft rows up/down in a block of CSC, used in the next func
 */
 void shift(int n_col, int n_row, double *a, int lda, const int inc);

/*
 * This function re-order the part of L facing to diagonal
 * according to given ordering and swap array.
 */
 void swapping_row_in_blockedL_new(int supWdt, int supNo, const size_t *lC, const size_t *Li_ptr, int *lR,
                                   const int curCol, const int *blockSet, double *lValues,
                                   int top, int *xi, int *lbs, int *ubs,
                                   int *cur_perm, int *cur_swap);


/*
 * Row reordering in blocked L matrix
 */
 void row_reordering(int supNo, const size_t *lC, const size_t *Li_ptr, int *lR,
                     const int *blockSet, int *aTree, int *cT, int *rT, int *col2Sup, double *lValues,
                     std::vector<int> swap_req, int *full_swap, int *xi, int *map, int *ws, double *wsf);

/*
 * The input is blocked diagonal matrix, the max block size is 2x2
 * Suitable for LDLT SBK .
 * Diagonal matrix is stored in array of size 2*n
 */

 void blocked_2by2_solver(int n, double *D, double *rhs,
                          int n_rhs, int lda, int lda_d);

/*
 * The input is blocked diagonal matrix, the max block size is 2x2
 * Suitable for LDLT SBK .
 * Diagonal matrix is stored in array of size 2*n
 */

 void blocked_2by2_solver_update(int n, double *D, double *rhs,
                                 int n_rhs, int lda, int lda_d, int *mask);

/*
 * The input is blocked diagonal matrix, the max block size is 2x2
 * Suitable for LDLT SBK .
 * Diagonal matrix is stored in array of size 2*n
 * diagonal matrix is symmetric
 */

 void blocked_2by2_mult(int n, int m, double *D, double *src, double *dst,
                        int lda, int lda_d);
}
#endif //PROJECT_BLASKERNELS_H
