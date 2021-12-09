//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/Etree.h"

#include "nasoq/common/SparseUtils.h"

namespace nasoq {

 void update_etree(int k, int i, int *Parent, int *Ancestor) {
  int a;
  for (;;)  /* traverse the path from k to the root of the tree */
  {
   a = Ancestor[k];
   if (a == i) {
    /* final ancestor reached; no change to tree */
    return;
   }
   /* perform path compression */
   Ancestor[k] = i;
   if (a == EMPTY) {
    /* final ancestor undefined; this is a new edge in the tree */
    Parent[k] = i;
    return;
   }
   /* traverse up to the ancestor of k */
   k = a;
  }
 }

 int etreeC(CSC *A, int *Parent, int &status) {
  int *Ap, *Ai, *Anz, *Ancestor, *Prev, *Iwork;
  int i, j, jprev, p, pend, nrow, ncol, packed, stype;
  size_t s;
  int ok = TRUE;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  /*RETURN_IF_NULL_COMMON (FALSE) ;
  RETURN_IF_NULL (A, FALSE) ;
  RETURN_IF_NULL (Parent, FALSE) ;
  RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;*/
  //Common->status = CHOLMOD_OK ;
  status = TRUE;

  /* ---------------------------------------------------------------------- */
  /* allocate workspace */
  /* ---------------------------------------------------------------------- */

  stype = A->stype;

  /* s = A->nrow + (stype ? 0 : A->ncol) */
  s = add_size_t(A->nrow, (stype ? 0 : A->ncol), &ok);
  if (!ok) {
   //ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
   return (FALSE);
  }

  //CHOLMOD(allocate_work) (0, s, 0, Common) ;
  /*if (Common->status < CHOLMOD_OK)
  {
      return (FALSE) ;	*//* out of memory *//*
    }*/

  //ASSERT (CHOLMOD(dump_sparse) (A, "etree", Common) >= 0) ;
  //Iwork = Common->Iwork ;
  Iwork = new int[s]();

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  ncol = A->ncol; /* the number of columns of A */
  nrow = A->nrow; /* the number of rows of A */
  Ap = A->p;  /* size ncol+1, column pointers for A */
  Ai = A->i;  /* the row indices of A */
  Anz = A->nz; /* number of nonzeros in each column of A */
  packed = A->packed;
  Ancestor = Iwork; /* size ncol (i/i/l) */

  for (j = 0; j < ncol; j++) {
   Parent[j] = EMPTY;
   Ancestor[j] = EMPTY;
  }

  /* ---------------------------------------------------------------------- */
  /* compute the etree */
  /* ---------------------------------------------------------------------- */

  if (stype > 0) {

   /* ------------------------------------------------------------------ */
   /* symmetric (upper) case: compute etree (A) */
   /* ------------------------------------------------------------------ */

   for (j = 0; j < ncol; j++) {
    /* for each row i in column j of triu(A), excluding the diagonal */
    p = Ap[j];
    pend = (packed) ? (Ap[j + 1]) : (p + Anz[j]);
    for (; p < pend; p++) {
     i = Ai[p];
     if (i < j) {
      update_etree(i, j, Parent, Ancestor);
     }
    }
   }

  } else if (stype == 0) {

   /* ------------------------------------------------------------------ */
   /* unsymmetric case: compute etree (A'*A) */
   /* ------------------------------------------------------------------ */

   Prev = Iwork + ncol; /* size nrow (i/i/l) */
   for (i = 0; i < nrow; i++) {
    Prev[i] = EMPTY;
   }
   for (j = 0; j < ncol; j++) {
    /* for each row i in column j of A */
    p = Ap[j];
    pend = (packed) ? (Ap[j + 1]) : (p + Anz[j]);
    for (; p < pend; p++) {
     /* a graph is constructed dynamically with one path per row
      * of A.  If the ith row of A contains column indices
      * (j1,j2,j3,j4) then the new graph has edges (j1,j2), (j2,j3),
      * and (j3,j4).  When at node i of this path-graph, all edges
      * (jprev,j) are considered, where jprev<j */
     i = Ai[p];
     jprev = Prev[i];
     if (jprev != EMPTY) {
      update_etree(jprev, j, Parent, Ancestor);
     }
     Prev[i] = j;
    }
   }

  } else {

   /* ------------------------------------------------------------------ */
   /* symmetric case with lower triangular part not supported */
   /* ------------------------------------------------------------------ */

   // ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
   return (FALSE);
  }

//    ASSERT (CHOLMOD(dump_parent) (Parent, ncol, "Parent", Common)) ;
  delete[]Iwork;
  return (TRUE);
 }

 int etree1(int ncol, int nrow, int *Ap, int *Ai, int *Anz, int stype, int packed, int *Parent) {
  int *Ancestor, *Prev, *Iwork;
  int i, j, jprev, p, pend;
  size_t s;
  int ok = TRUE;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

/*    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;*/

  /* ---------------------------------------------------------------------- */
  /* allocate workspace */
  /* ---------------------------------------------------------------------- */

  //stype = A->stype ;

  /* s = A->nrow + (stype ? 0 : A->ncol) */
/*    s = CHOLMOD(add_size_t) (A->nrow, (stype ? 0 : A->ncol), &ok) ;
    if (!ok)
    {
        ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
        return (FALSE) ;
    }

    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
        return (FALSE) ;	*//* out of memory *//*
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "etree", Common) >= 0) ;
    Iwork = Common->Iwork ;*/

  s = nrow + (stype ? 0 : ncol);
  Iwork = new int[s]();
  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  //ncol = A->ncol ;	/* the number of columns of A */
  //nrow = A->nrow ;	/* the number of rows of A */
  //Ap = A->p ;		/* size ncol+1, column pointers for A */
  //Ai = A->i ;		/* the row indices of A */
  //Anz = A->nz ;	/* number of nonzeros in each column of A */
  //packed = A->packed ;
  Ancestor = Iwork; /* size ncol (i/i/l) */

  for (j = 0; j < ncol; j++) {
   Parent[j] = EMPTY;
   Ancestor[j] = EMPTY;
  }

  /* ---------------------------------------------------------------------- */
  /* compute the etree */
  /* ---------------------------------------------------------------------- */

  if (stype > 0) {

   /* ------------------------------------------------------------------ */
   /* symmetric (upper) case: compute etree (A) */
   /* ------------------------------------------------------------------ */

   for (j = 0; j < ncol; j++) {
    /* for each row i in column j of triu(A), excluding the diagonal */
    p = Ap[j];
    pend = (packed) ? (Ap[j + 1]) : (p + Anz[j]);
    for (; p < pend; p++) {
     i = Ai[p];
     if (i < j) {
      update_etree(i, j, Parent, Ancestor);
     }
    }
   }

  } else if (stype == 0) {

   /* ------------------------------------------------------------------ */
   /* unsymmetric case: compute etree (A'*A) */
   /* ------------------------------------------------------------------ */

   Prev = Iwork + ncol; /* size nrow (i/i/l) */
   for (i = 0; i < nrow; i++) {
    Prev[i] = EMPTY;
   }
   for (j = 0; j < ncol; j++) {
    /* for each row i in column j of A */
    p = Ap[j];
    pend = (packed) ? (Ap[j + 1]) : (p + Anz[j]);
    for (; p < pend; p++) {
     /* a graph is constructed dynamically with one path per row
      * of A.  If the ith row of A contains column indices
      * (j1,j2,j3,j4) then the new graph has edges (j1,j2), (j2,j3),
      * and (j3,j4).  When at node i of this path-graph, all edges
      * (jprev,j) are considered, where jprev<j */
     i = Ai[p];
     jprev = Prev[i];
     if (jprev != EMPTY) {
      update_etree(jprev, j, Parent, Ancestor);
     }
     Prev[i] = j;
    }
   }

  } else {

   /* ------------------------------------------------------------------ */
   /* symmetric case with lower triangular part not supported */
   /* ------------------------------------------------------------------ */

   printf("symmetric lower not supported");
   return (FALSE);
  }

  //ASSERT (CHOLMOD(dump_parent) (Parent, ncol, "Parent", Common)) ;
  return (TRUE);
 }

 int *etree(int n, int *Ap, int *Ai, int ata) {
  /* compute the etree of A (using triu(A),
   * or A'A without forming A'A
   * n: the matrix size
   * Ap: column pointer
   * Ai: row index
   * ata: A'A or A
   * */
  int i, k, p, m, inext, *w, *parent, *ancestor, *prev;
  if (n < 0 || Ap == NULL || Ai == NULL) //check inputs
   return 0;
  m = n;
  parent = new int[n]; //result allocation
  w = new int[n + (ata ? m : 0)]; // get workspace
  if (w == NULL || parent == NULL)
   return 0;

  ancestor = w;
  prev = w + n;
  if (ata) for (i = 0; i < m; i++) prev[i] = -1;
  for (k = 0; k < n; k++) {
   parent[k] = -1;                   /* node k has no parent yet */
   ancestor[k] = -1;                 /* nor does k have an ancestor */
   for (p = Ap[k]; p < Ap[k + 1]; p++) {
    i = ata ? (prev[Ai[p]]) : (Ai[p]);
    for (; i != -1 && i < k; i = inext)   /* traverse from i to k */
    {
     inext = ancestor[i];              /* inext = ancestor of i */
     ancestor[i] = k;                  /* path compression */
     if (inext == -1) parent[i] = k;   /* no anc., parent is k */
    }
    if (ata) prev[Ai[p]] = k;
   }
  }
  delete[]w;
  return parent;
 }
}