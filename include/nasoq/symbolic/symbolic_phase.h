//
// Created by kazem on 4/11/19.
//

#ifndef PROJECT_SYMBOLIC_PHASE_H
#define PROJECT_SYMBOLIC_PHASE_H

#include <cstddef>

#include "nasoq/common/def.h"

#ifdef SCOTCH
#include "scotch.h"
#endif
#ifdef SUITESPARSE_FOUND
#include "amd.h"
#endif
#include "metis.h"
#include "InspectionLevel_06.h"
#include "performanceModel.h"

namespace nasoq {
//#undef METIS
 int permute_matrices
   (
     /* ---- input ---- */
     CSC *A, /* matrix to permute */
     int ordering, /* ordering method used */
     int *Perm,  /* fill-reducing permutation */
     int *fset,  /* subset of 0:(A->ncol)-1 */
     size_t fsize, /* size of fset */
     int do_rowcolcounts,/* if TRUE, compute both S and F.  If FALSE, only
			 * S is needed for the symmetric case, and only F for
			 * the unsymmetric case */

     /* ---- output --- */
     CSC **A1_handle,     /* see comments below for A1, A2, S, F */
     CSC **A2_handle,
     CSC **S_handle,
     CSC **F_handle,
     /* --------------- */
     //cholmod_common *Common
     int &status
   );

/*
 *
 */

 int analyze_ordering
   (
     /* ---- input ---- */
     CSC *A, /* matrix to analyze */
     int ordering, /* ordering method used */
     int *Perm,  /* size n, fill-reducing permutation to analyze */
     int *fset,  /* subset of 0:(A->ncol)-1 */
     size_t fsize, /* size of fset */
     /* ---- output --- */
     int *Parent, /* size n, elimination tree */
     int *Post,  /* size n, postordering of elimination tree */
     int *ColCount, /* size n, nnz in each column of L */

     /* ---- workspace  */
     int *First,  /* size n workspace for cholmod_postorder */
     int *Level,  /* size n workspace for cholmod_postorder */
     /* --------------- */
     // cholmod_common *Common
     int &status
   );


/* Ordering and analysis for sparse Cholesky or sparse QR.
 * */

 BCSC *symbolic_analysis_lin_solve
   (
     /* ---- input ---- */
     int for_whom,       /* FOR_SPQR     (0): for SPQR but not GPU-accelerated
                           FOR_CHOLESKY (1): for Cholesky (GPU or not)
                           FOR_SPQRGPU  (2): for SPQR with GPU acceleration */
     CSC *A, /* matrix to order and analyze */
     int *UserPerm, /* user-provided permutation, size A->nrow */
     int *fset,  /* subset of 0:(A->ncol)-1 */
     int *nrelax,
     double *zrelax,
     size_t fsize, /* size of fset */
     int *&prunePtr,
     int *&pruneSet,
     int &levelNo,
     int *&levelPtr,
     int *&levelSet,
     int &parNo,
     int *&parPtr,
     int *&partition,
     int &levelNo_s,
     int *&levelPtr_s,
     int &parNo_s,
     int *&parPtr_s,
     int *&partition_s,
     int costParam,
     int levelParam, //min level distance
     int finalSeqNodes, //number of level in final partition
     /* --------------- */
     int status,
     //out
     int &maxSupWid,
     int &maxCol,
     double &orderingTime,
     int simplicial_en,
     int *extra_cols,
     size_t *inPerm = NULL
   );

}

#endif //PROJECT_SYMBOLIC_PHASE_H
