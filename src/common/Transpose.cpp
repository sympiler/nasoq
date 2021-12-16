//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/Transpose.h"

#include <cstdlib>

namespace nasoq {

 int transpose_sym_real(CSC *A, int *Perm, CSC *F, int *Wi, int *Pinv, int &status) {
  double *Ax, *Az, *Fx, *Fz;
  int *Ap, *Anz, *Ai, *Fp, *Fj, *Iwork;
  int p, pend, packed, fp, upper, permute, jold, n, i, j, iold;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  /* ensure the xtype of A and F match (ignored if this is pattern version) */
/*    if (!XTYPE_OK (A->xtype))
    {
       // ERROR (CHOLMOD_INVALID, "real/complex mismatch") ;
        return (FALSE) ;
    }*/

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  permute = (Perm != NULL);
  n = A->nrow;
  Ap = A->p;  /* size A->ncol+1, column pointers of A */
  Ai = A->i;  /* size nz = Ap [A->ncol], row indices of A */
  Ax = A->x;  /* size nz, real values of A */
  Az = A->z;  /* size nz, imag values of A */
  Anz = A->nz;
  packed = A->packed;
  // ASSERT (IMPLIES (!packed, Anz != NULL)) ;
  upper = (A->stype > 0);

  Fp = F->p;  /* size A->nrow+1, row pointers of F */
  Fj = F->i;  /* size nz, column indices of F */
  Fx = F->x;  /* size nz, real values of F */
  Fz = F->z;  /* size nz, imag values of F */

  /* ---------------------------------------------------------------------- */
  /* get workspace */
  /* ---------------------------------------------------------------------- */

  //Iwork = Common->Iwork ;
  // Iwork = new int[2*n]();
  // Wi = Iwork ;	/* size n (i/l/l) */
  // Pinv = Iwork + n ;	/* size n (i/i/l) , unused if Perm NULL */

  /* ---------------------------------------------------------------------- */
  /* construct the transpose */
  /* ---------------------------------------------------------------------- */

  if (permute) {
   if (upper) {
    /* permuted, upper */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     p = Ap[jold];
     pend = (packed) ? Ap[jold + 1] : p + Anz[jold];
     for (; p < pend; p++) {
      iold = Ai[p];
      if (iold <= jold) {
       i = Pinv[iold];
       if (i < j) {
        fp = Wi[i]++;
        Fj[fp] = j;
        Fx[fp] = Ax[p];
       } else {
        fp = Wi[j]++;
        Fj[fp] = i;
        Fx[fp] = Ax[p];
       }
      }
     }
    }
   } else {
    /* permuted, lower */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     p = Ap[jold];
     pend = (packed) ? Ap[jold + 1] : p + Anz[jold];
     for (; p < pend; p++) {
      iold = Ai[p];
      if (iold >= jold) {
       i = Pinv[iold];
       if (i > j) {
        fp = Wi[i]++;
        Fj[fp] = j;
        Fx[fp] = Ax[p];
       } else {
        fp = Wi[j]++;
        Fj[fp] = i;
        Fx[fp] = Ax[p];
       }
      }
     }
    }
   }
  } else {
   if (upper) {
    /* unpermuted, upper */
    for (j = 0; j < n; j++) {
     p = Ap[j];
     pend = (packed) ? Ap[j + 1] : p + Anz[j];
     for (; p < pend; p++) {
      i = Ai[p];
      if (i <= j) {
       fp = Wi[i]++;
       Fj[fp] = j;
       Fx[fp] = Ax[p];
      }
     }
    }
   } else {
    /* unpermuted, lower */
    for (j = 0; j < n; j++) {
     p = Ap[j];
     pend = (packed) ? Ap[j + 1] : p + Anz[j];
     for (; p < pend; p++) {
      i = Ai[p];
      if (i >= j) {
       fp = Wi[i]++;
       Fj[fp] = j;
       Fx[fp] = Ax[p];
      }
     }
    }
   }
  }

  return (TRUE);
 }

 int transpose_sym(CSC *A, int values, int *Perm, CSC *F, int status) {
  int *Ap, *Anz, *Ai, *Fp, *Wi, *Pinv, *Iwork;
  int p, pend, packed, upper, permute, jold, n, i, j, k, iold;
  size_t s;
  int ok = TRUE;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

/*    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (F, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;*/
  if (A->nrow != A->ncol || A->stype == 0) {
   /* this routine handles square symmetric matrices only */
   //ERROR (CHOLMOD_INVALID, "matrix must be symmetric") ;
   return (FALSE);
  }
  if (A->nrow != F->ncol || A->ncol != F->nrow) {
   //    ERROR (CHOLMOD_INVALID, "F has the wrong dimensions") ;
   return (FALSE);
  }
  // Common->status = CHOLMOD_OK ;
  status = TRUE;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  permute = (Perm != NULL);
  n = A->nrow;
  Ap = A->p;  /* size A->ncol+1, column pointers of A */
  Ai = A->i;  /* size nz = Ap [A->ncol], row indices of A */
  Anz = A->nz;
  packed = A->packed;
  // ASSERT (IMPLIES (!packed, Anz != NULL)) ;
  upper = (A->stype > 0);

  Fp = F->p;  /* size A->nrow+1, row pointers of F */

  /* ---------------------------------------------------------------------- */
  /* allocate workspace */
  /* ---------------------------------------------------------------------- */

  /* s = (Perm != NULL) ? 2*n : n */
  s = add_size_t(n, ((Perm != NULL) ? n : 0), &ok);
  if (!ok) {
   //    ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
   return (FALSE);
  }

  /*CHOLMOD(allocate_work) (0, s, 0, Common) ;
  if (Common->status < CHOLMOD_OK)
  {
      return (FALSE) ;	*//* out of memory *//*
    }*/

  /* ---------------------------------------------------------------------- */
  /* get workspace */
  /* ---------------------------------------------------------------------- */

  //Iwork = Common->Iwork ;
  Iwork = (int *) calloc(s, sizeof(int));//new int[s]();
  Wi = Iwork;     /* size n (i/l/l) */
  Pinv = Iwork + n;     /* size n (i/i/l) , unused if Perm NULL */

  /* ---------------------------------------------------------------------- */
  /* check Perm and construct inverse permutation */
  /* ---------------------------------------------------------------------- */

  if (permute) {
   for (i = 0; i < n; i++) {
    Pinv[i] = EMPTY;
   }
   for (k = 0; k < n; k++) {
    i = Perm[k];
    if (i < 0 || i > n || Pinv[i] != EMPTY) {
     //        ERROR (CHOLMOD_INVALID, "invalid permutation") ;
     return (FALSE);
    }
    Pinv[i] = k;
   }
  }

  /* Perm is now valid */
  //ASSERT (CHOLMOD(dump_perm) (Perm, n, n, "Perm", Common)) ;

  /* ---------------------------------------------------------------------- */
  /* count the entries in each row of F */
  /* ---------------------------------------------------------------------- */

  for (i = 0; i < n; i++) {
   Wi[i] = 0;
  }

  if (packed) {
   if (permute) {
    if (upper) {
     /* packed, permuted, upper */
     for (j = 0; j < n; j++) {
      jold = Perm[j];
      pend = Ap[jold + 1];
      for (p = Ap[jold]; p < pend; p++) {
       iold = Ai[p];
       if (iold <= jold) {
        i = Pinv[iold];
        Wi[MIN (i, j)]++;
       }
      }
     }
    } else {
     /* packed, permuted, lower */
     for (j = 0; j < n; j++) {
      jold = Perm[j];
      pend = Ap[jold + 1];
      for (p = Ap[jold]; p < pend; p++) {
       iold = Ai[p];
       if (iold >= jold) {
        i = Pinv[iold];
        Wi[MAX (i, j)]++;
       }
      }
     }
    }
   } else {
    if (upper) {
     /* packed, unpermuted, upper */
     for (j = 0; j < n; j++) {
      pend = Ap[j + 1];
      for (p = Ap[j]; p < pend; p++) {
       i = Ai[p];
       if (i <= j) {
        Wi[i]++;
       }
      }
     }
    } else {
     /* packed, unpermuted, lower */
     for (j = 0; j < n; j++) {
      pend = Ap[j + 1];
      for (p = Ap[j]; p < pend; p++) {
       i = Ai[p];
       if (i >= j) {
        Wi[i]++;
       }
      }
     }
    }
   }
  } else {
   if (permute) {
    if (upper) {
     /* unpacked, permuted, upper */
     for (j = 0; j < n; j++) {
      jold = Perm[j];
      p = Ap[jold];
      pend = p + Anz[jold];
      for (; p < pend; p++) {
       iold = Ai[p];
       if (iold <= jold) {
        i = Pinv[iold];
        Wi[MIN (i, j)]++;
       }
      }
     }
    } else {
     /* unpacked, permuted, lower */
     for (j = 0; j < n; j++) {
      jold = Perm[j];
      p = Ap[jold];
      pend = p + Anz[jold];
      for (; p < pend; p++) {
       iold = Ai[p];
       if (iold >= jold) {
        i = Pinv[iold];
        Wi[MAX (i, j)]++;
       }
      }
     }
    }
   } else {
    if (upper) {
     /* unpacked, unpermuted, upper */
     for (j = 0; j < n; j++) {
      p = Ap[j];
      pend = p + Anz[j];
      for (; p < pend; p++) {
       i = Ai[p];
       if (i <= j) {
        Wi[i]++;
       }
      }
     }
    } else {
     /* unpacked, unpermuted, lower */
     for (j = 0; j < n; j++) {
      p = Ap[j];
      pend = p + Anz[j];
      for (; p < pend; p++) {
       i = Ai[p];
       if (i >= j) {
        Wi[i]++;
       }
      }
     }
    }
   }
  }

  /* ---------------------------------------------------------------------- */
  /* compute the row pointers */
  /* ---------------------------------------------------------------------- */

  p = 0;
  for (i = 0; i < n; i++) {
   Fp[i] = p;
   p += Wi[i];
  }
  Fp[n] = p;
  for (i = 0; i < n; i++) {
   Wi[i] = Fp[i];
  }

  if (p > (int) (F->nzmax)) {
   //  ERROR (CHOLMOD_INVALID, "F is too small") ;
   return (FALSE);
  }

  /* ---------------------------------------------------------------------- */
  /* transpose matrix, using template routine */
  /* ---------------------------------------------------------------------- */

  ok = FALSE;
  if (values == 0 || F->xtype == CHOLMOD_PATTERN) {
   //  PRINT2 (("\n:::: p_transpose_sym Perm %p\n", Perm)) ;
   ok = transpose_sym_real(A, Perm, F, Wi, Pinv, status);
  } else if (F->xtype == CHOLMOD_REAL) {
   // PRINT2 (("\n:::: r_transpose_sym Perm %p\n", Perm)) ;
   ok = transpose_sym_real(A, Perm, F, Wi, Pinv, status);
  }
///
  /* ---------------------------------------------------------------------- */
  /* finalize result F */
  /* ---------------------------------------------------------------------- */

  /* F is sorted if there is no permutation vector */
  if (ok) {
   F->sorted = !permute;
   F->packed = TRUE;
   F->stype = -SIGN (A->stype); /* flip the stype */
   //  ASSERT (CHOLMOD(dump_sparse) (F, "output F sym", Common) >= 0) ;
  }
  free(Iwork);
  return (ok);
 }

 CSC *ptranspose(CSC *A, int values, int *Perm, int *fset, size_t fsize, int status) {
  int *Ap, *Anz;
  CSC *F;
  int nrow, ncol, use_fset, j, jj, fnz, packed, stype, nf, xtype;
  size_t ineed;
  int ok = TRUE;

  nf = fsize;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  /*RETURN_IF_NULL_COMMON (NULL) ;
  RETURN_IF_NULL (A, FALSE) ;
  RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;*/
  stype = A->stype;
  //Common->status = CHOLMOD_OK ;
  status = TRUE;

  /* ---------------------------------------------------------------------- */
  /* allocate workspace */
  /* ---------------------------------------------------------------------- */

  nrow = A->nrow;
  ncol = A->ncol;

  if (stype != 0) {
   use_fset = FALSE;
   if (Perm != NULL) {
    ineed = mult_size_t(A->nrow, 2, &ok);
   } else {
    ineed = A->nrow;
   }
  } else {
   use_fset = (fset != NULL);
   if (use_fset) {
    ineed = MAX (A->nrow, A->ncol);
   } else {
    ineed = A->nrow;
   }
  }

  if (!ok) {
   //ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
   return (NULL);
  }

  /*CHOLMOD(allocate_work) (0, ineed, 0, Common) ;
  if (Common->status < CHOLMOD_OK)
  {
      return (NULL) ;	    *//* out of memory *//*
    }*/

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  Ap = A->p;
  Anz = A->nz;
  packed = A->packed;
  // ASSERT (IMPLIES (!packed, Anz != NULL)) ;
  xtype = values ? A->xtype : CHOLMOD_PATTERN;

  /* ---------------------------------------------------------------------- */
  /* allocate F */
  /* ---------------------------------------------------------------------- */

  /* determine # of nonzeros in F */
  if (stype != 0) {
   /* F=A' or F=A(p,p)', fset is ignored */
   fnz = getNNZ(A->ncol, A->p, A->nz, A->packed, status);
  } else {
   nf = (use_fset) ? nf : ncol;
   if (use_fset) {
    fnz = 0;
    /* F=A(:,f)' or F=A(p,f)' */
    for (jj = 0; jj < nf; jj++) {
     /* The fset is not yet checked; it will be thoroughly checked
      * in cholmod_transpose_unsym.  For now, just make sure we don't
      * access Ap and Anz out of bounds. */
     j = fset[jj];
     if (j >= 0 && j < ncol) {
      fnz += packed ? (Ap[j + 1] - Ap[j]) : MAX (0, Anz[j]);
     }
    }
   } else {
    /* F=A' or F=A(p,:)' */
    //fnz = CHOLMOD(nnz) (A, Common) ;
    fnz = getNNZ(A->ncol, A->p, A->nz, A->packed, status);
   }
  }

  /* F is ncol-by-nrow, fnz nonzeros, sorted unless f is present and unsorted,
   * packed, of opposite stype as A, and with/without numerical values */
  /* F = CHOLMOD(allocate_sparse) (ncol, nrow, fnz, TRUE, TRUE, -SIGN(stype),
                                 xtype, Common) ;*/ //TODO
  F = new CSC;
  //F = (CSC*) calloc(1 , sizeof(CSC));
  allocateAC(F, A->nrow, fnz, -SIGN(A->stype), TRUE);

  if (!status) {
   return (NULL);     /* out of memory */
  }

  /* ---------------------------------------------------------------------- */
  /* transpose and optionally permute the matrix A */
  /* ---------------------------------------------------------------------- */

  if (stype != 0) {
   /* F = A (p,p)', using upper or lower triangular part of A only */
   ok = transpose_sym(A, values, Perm, F, status);
  } else {
   /* F = A (p,f)' */
   // ok = CHOLMOD(transpose_unsym) (A, values, Perm, fset, nf, F, Common) ;
   printf("unsym is not supported.");
  }

  /* ---------------------------------------------------------------------- */
  /* return the matrix F, or NULL if an error occured */
  /* ---------------------------------------------------------------------- */

  /*  if (!ok)
    {
        CHOLMOD(free_sparse) (&F, Common) ;
    }*/
  return (F);
 }

 CSC *transposeC(CSC *A, int values, int &status) {
  return (ptranspose(A, values, NULL, NULL, 0, status));
 }

 double cumsum(int *p, int *c, int n) {
  int i, nz = 0;
  double nz2 = 0;
  if (!p || !c) return (-1);     /* check inputs */
  for (i = 0; i < n; i++) {
   p[i] = nz;
   nz += c[i];
   nz2 += c[i];              /* also in double to avoid csi overflow */
   c[i] = p[i];             /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p[n] = nz;
  return (nz2);                  /* return sum (c [0..n-1]) */
 }

 int transpose(int n, int m, int *Ap, int *Ai, double *Ax, int values, int *Cp, int *Ci, double *Cx) {
  int p, q, j, *w;
  //double *Cx;
  //cs *C ;
  if (!Ai || !Ap) return (0);    /* check inputs */
  //m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  //C = cs_spalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
  //Cp = new int[n+1]; Ci = new int[Ap[n]];
  //if(values) Cx = new double[Ap[n]];
  //else Cx=NULL;
  //w = cs_calloc (m, sizeof (csi)) ;                      /* get workspace */
  w = new int[m]();
  if (!Cp || !w) return -1;       /* out of memory */
  //Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (p = 0; p < Ap[n]; p++) w[Ai[p]]++;          /* row counts */
  cumsum(Cp, w, m);                                 /* row pointers */
  for (j = 0; j < n; j++) {
   for (p = Ap[j]; p < Ap[j + 1]; p++) {
    Ci[q = w[Ai[p]]++] = j; /* place A(i,j) as entry C(j,i) */
    if (Cx) Cx[q] = Ax[p];
   }
  }
  delete[]w;
  return 1;  /* success;  */
 }
}