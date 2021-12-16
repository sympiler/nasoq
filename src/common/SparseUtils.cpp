//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/SparseUtils.h"

#include <cmath>

namespace nasoq {

 long int getNNZ(int ncol, int *Ap, int *Anz, int packed, int &status) {
  size_t nz;
  int j;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  /*RETURN_IF_NULL_COMMON (EMPTY) ;
  RETURN_IF_NULL (A, EMPTY) ;
  RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
  Common->status = CHOLMOD_OK ;*/
  status = TRUE;
  /* ---------------------------------------------------------------------- */
  /* return nnz (A) */
  /* ---------------------------------------------------------------------- */

  if (packed) {
   //Ap = A->p ;
   //RETURN_IF_NULL (Ap, EMPTY) ;
   nz = Ap[ncol];
  } else {
   //Anz = A->nz ;
   //RETURN_IF_NULL (Anz, EMPTY) ;
   nz = 0;
   for (j = 0; j < ncol; j++) {
    nz += MAX (0, Anz[j]);
   }
  }
  return (nz);
 }

 size_t add_size_t(size_t a, size_t b, int *ok) {
  size_t s = a + b;
  (*ok) = (*ok) && (s >= a);
  return ((*ok) ? s : 0);
 }

 size_t mult_size_t(size_t a, size_t k, int *ok) {
  size_t p = 0, s;
  while (*ok) {
   if (k % 2) {
    p = p + a;
    (*ok) = (*ok) && (p >= a);
   }
   k = k / 2;
   if (!k) return (p);
   s = a + a;
   (*ok) = (*ok) && (s >= a);
   a = s;
  }
  return (0);
 }

 int makeUnique(int *node2Par, std::vector<int> &list, int n, bool *ws) {
  int min = INT32_MAX;
  for (int ii = 0; ii < list.size();) {
   int tmp = node2Par[list[ii]];
   if (!ws[tmp]) {//if first time
    ws[tmp] = true;
    min = min > tmp ? tmp : min;
    ii++;
   } else {//otherwise remove it
    list.erase(list.begin() + ii);
   }
  }
//Reset it for future use
  for (int ii = 0; ii < list.size(); ++ii) {
   int tmp = node2Par[list[ii]];
   ws[tmp] = false;
  }
  return min;//returns cluster with min number.
 }

 void add_vec(int n, double *a, double alpha, double *b) {
  //std::fill(b, b+n, 0.0);
  if (alpha == 0.0)
   return;
  if (alpha == 1.0)
   for (int i = 0; i < n; ++i) {
    double tmp = *b + *(a++);
    *(b++) = tmp;
   }
  else if (alpha == -1.0) {
   for (int i = 0; i < n; ++i) {
    double tmp = *b - *(a++);
    *(b++) = tmp;
   }
  } else
   for (int i = 0; i < n; ++i) {
    double tmp = *b + alpha * *(a++);
    *(b++) = tmp;
   }
 }

 void sparse2dense(CSC *A, double *D) {
  std::fill_n(D, A->ncol * A->nrow, 0);
  for (int i = 0; i < A->ncol; ++i) {
   for (int j = A->p[i]; j < A->p[i + 1]; ++j) {
    int r = A->i[j];
    D[i * A->nrow + r] = A->x[j];
    if (A->stype == -1 || A->stype == 1) {
     D[r * A->ncol + i] = A->x[j];
    }
   }
  }
 }

 void dense_traspose(int r, int c, double *a_in, double *a_out) {
  for (int i = 0; i < c; ++i) {
   for (int j = 0; j < r; ++j) {
    a_out[i * r + j] = a_in[j * c + i]; //a[i,j] = a[j,i];
   }
  }
 }

 int compute_recieporical_length(CSC *A, double *a, double *rec_norm) {
  int stype = A->stype;
  int *Ap = A->p;
  int ncol = A->ncol;
  int *Ai = A->i;
  double *Ax = A->x;
  std::fill_n(rec_norm, A->nrow, 0.0);
  if (stype == 0) {
   // compute square of rows
   for (int j = 0; j < ncol; j++) {
    for (int k = Ap[j]; k < Ap[j + 1]; ++k) {
     int p = Ai[k];
     double s = Ax[k];
     rec_norm[p] += (s * s);
    }
   }
/*  for (int l = 0; l < ncol; ++l) {
   std::cout<<rec_norm[l]<<"\n";
  }*/
   for (int i = 0; i < A->nrow; ++i) {
    if (rec_norm[i] > 0) {
     rec_norm[i] = 1.0 / std::sqrt(rec_norm[i]);
    } else {
     if (a[i] != 0)
      return 0;
    }
   }
  }
/* for (int l = 0; l < ncol; ++l) {
  std::cout<<rec_norm[l]<<"\n";
 }*/
  return 1;
 }
}