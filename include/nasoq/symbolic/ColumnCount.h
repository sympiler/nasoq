//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_COLUMNCOUNT_H
#define CHOLOPENMP_COLUMNCOUNT_H

#include <cstddef>

#include "nasoq/common/def.h"

namespace nasoq {

/* ========================================================================== */
/* === initialize_node ====================================================== */
/* ========================================================================== */

 int initialize_node  /* initial work for kth node in postordered etree */
   (
     int k,  /* at the kth step of the algorithm (and kth node) */
     int Post[], /* Post [k] = i, the kth node in postordered etree */
     int Parent[], /* Parent [i] is the parent of i in the etree */
     int ColCount[], /* ColCount [c] is the current weight of node c */
     int PrevNbr[] /* PrevNbr [u] = k if u was last considered at step k */
   );


/* ========================================================================== */
/* === process_edge ========================================================= */
/* ========================================================================== */

/* edge (p,u) is being processed.  p < u is a descendant of its ancestor u in
 * the etree.  node p is the kth node in the postordered etree.  */

 void process_edge
   (
     int p,  /* process edge (p,u) of the matrix */
     int u,
     int k,  /* we are at the kth node in the postordered etree */
     int First[], /* First [i] = k if the postordering of first
			 * descendent of node i is k */
     int PrevNbr[], /* u was last considered at step k = PrevNbr [u] */
     int ColCount[], /* ColCount [c] is the current weight of node c */
     int PrevLeaf[], /* s = PrevLeaf [u] means that s was the last leaf
			 * seen in the subtree rooted at u.  */
     int RowCount[], /* RowCount [i] is # of nonzeros in row i of L,
			 * including the diagonal.  Not computed if NULL. */
     int SetParent[], /* the FIND/UNION data structure, which forms a set
			 * of trees.  A root i has i = SetParent [i].  Following
			 * a path from i to the root q of the subtree containing
			 * i means that q is the SetParent representative of i.
			 * All nodes in the tree could have their SetParent
			 * equal to the root q; the tree representation is used
			 * to save time.  When a path is traced from i to its
			 * root q, the path is re-traversed to set the SetParent
			 * of the whole path to be the root q. */
     int Level[]  /* Level [i] = length of path from node i to root */
   );


/* ========================================================================== */
/* === finalize_node ======================================================== */
/* ========================================================================== */

 void finalize_node    /* compute UNION (p, Parent [p]) */
   (
     int p,
     int Parent[], /* Parent [p] is the parent of p in the etree */
     int SetParent[] /* see process_edge, above */
   );

/* ========================================================================== */
/* === cholmod_rowcolcounts ================================================= */
/* ========================================================================== */

 int rowcolcounts(
   /* ---- input ---- */
   CSC *A, /* matrix to analyze */
   int *fset,  /* subset of 0:(A->ncol)-1 */
   size_t fsize, /* size of fset */
   int *Parent, /* size nrow.  Parent [i] = p if p is the parent of i */
   int *Post,  /* size nrow.  Post [k] = i if i is the kth node in
			 * the postordered etree. */
   /* ---- output --- */
   int *RowCount, /* size nrow. RowCount [i] = # entries in the ith row of
			 * L, including the diagonal. */
   int *ColCount, /* size nrow. ColCount [i] = # entries in the ith
			 * column of L, including the diagonal. */
   int *First,  /* size nrow.  First [i] = k is the least postordering
			 * of any descendant of i. */
   int *Level,  /* size nrow.  Level [i] is the length of the path from
			 * i to the root, with Level [root] = 0. */
   /* --------------- */
   int &fl,
   int &aatfl,
   int &lnz,
   //   cholmod_common *Common
   int status
 );

/* ========================================================================== */
/* === cholmod_rowcolcounts ================================================= */
/* ========================================================================== */

 int rowcolcounts1
   (
     /* ---- input ---- */
     //cholmod_sparse *A,	/* matrix to analyze */
     int Anrow,
     int Ancol,
     int *Ap,
     int *Ai,
     int *Anz,
     double *Ax,
     int Astype,
     int packed,
     int *fset,  /* subset of 0:(A->ncol)-1 */
     size_t fsize, /* size of fset */
     int *Parent, /* size nrow.  Parent [i] = p if p is the parent of i */
     int *Post,  /* size nrow.  Post [k] = i if i is the kth node in
			 * the postordered etree. */
     /* ---- output --- */
     int *RowCount, /* size nrow. RowCount [i] = # entries in the ith row of
			 * L, including the diagonal. */
     int *ColCount, /* size nrow. ColCount [i] = # entries in the ith
			 * column of L, including the diagonal. */
     int *First,  /* size nrow.  First [i] = k is the least postordering
			 * of any descendant of i. */
     int *Level,  /* size nrow.  Level [i] is the length of the path from
			 * i to the root, with Level [root] = 0. */
     /* --------------- */
     int &aatfl,
     int &lnz,
     int &fl,
     // cholmod_common *Common
     int &status
   );
}
#endif //CHOLOPENMP_COLUMNCOUNT_H
