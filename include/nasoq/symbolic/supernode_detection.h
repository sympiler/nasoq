//
// Created by kazem on 5/14/19.
//

#ifndef PARS_SUPERNODE_DETECTION_H
#define PARS_SUPERNODE_DETECTION_H
/* ========================================================================== */
/* === subtree ============================================================== */
/* ========================================================================== */

#include "nasoq/common/def.h"

namespace nasoq {

/* In the symmetric case, traverse the kth row subtree from the nonzeros in
 * A (0:k1-1,k) and add the new entries found to the pattern of the kth row
 * of L.  The current supernode s contains the diagonal block k1:k2-1, so it
 * can be skipped.
 *
 * In the unsymmetric case, the nonzero pattern of A*F is computed one column
 * at a time (thus, the total time spent in this function is bounded below by
 * the time taken to multiply A*F, which can be high if A is tall and thin).
 * The kth column is A*F(:,k), or the set union of all columns A(:,j) for which
 * F(j,k) is nonzero.  This routine is called once for each entry j.  Only the
 * upper triangular part is needed, so only A (0:k1-1,j) is accessed, where
 * k1:k2-1 are the columns of the current supernode s (k is in the range k1 to
 * k2-1).
 *
 * If A is sorted, then the total time taken by this function is proportional
 * to the number of nonzeros in the strictly block upper triangular part of A,
 * plus the number of entries in the strictly block lower triangular part of
 * the supernodal part of L.  This excludes entries in the diagonal blocks
 * corresponding to the columns in each supernode.  That is, if k1:k2-1 are
 * in a single supernode, then only A (0:k1-1,k1:k2-1) are accessed.
 *
 * For the unsymmetric case, only the strictly block upper triangular part
 * of A*F is constructed.
 *
 * Only adds column indices corresponding to the leading columns of each
 * relaxed supernode.
 */

 void subtree
   (
     /* inputs, not modified: */
     int j,  /* j = k for symmetric case */
     int k,
     int Ap[],
     int Ai[],
     int Anz[],
     int SuperMap[],
     int Sparent[],
     int mark,
     int sorted,         /* true if the columns of A are sorted */
     int k1,             /* only consider A (0:k1-1,k) */

     /* input/output: */
     int Flag[],
     int Ls[],
     int Lpi2[]
   );


/* clear workspace used by cholmod_super_symbolic */
#define FREE_WORKSPACE \
{ \
    /* CHOLMOD(clear_flag) (Common) ; */ \
    CHOLMOD_CLEAR_FLAG (Common) ; \
    for (k = 0 ; k <= nfsuper ; k++) \
    { \
 Head [k] = EMPTY ; \
    } \
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ; \
} \


/* ========================================================================== */
/* === cholmod_super_symbolic2 ============================================== */
/* ========================================================================== */

/* Analyze for supernodal Cholesky or multifrontal QR. */

 int super_node_detection
   (
     /* ---- input ---- */
     int for_whom,       /* FOR_SPQR     (0): for SPQR but not GPU-accelerated
                           FOR_CHOLESKY (1): for Cholesky (GPU or not)
                           FOR_SPQRGPU  (2): for SPQR with GPU acceleration */
     CSC *A, /* matrix to analyze */
     CSC *F, /* F = A' or A(:,f)' */
     int *Parent, /* elimination tree */
     int *orig_or_extra, /*if 1 the col is update/downdate*/
     /* ---- in/out --- */
     BCSC *L, /* simplicial symbolic on input,
			 * supernodal symbolic on output */
     /* --------------- */
     int *nrelax,
     double *zrelax,
     // cholmod_common *Common
     int *Sparent,
     int status
   );
}
#endif //PARS_SUPERNODE_DETECTION_H
