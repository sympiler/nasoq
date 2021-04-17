//
// Created by kazem on 10/9/19.
//

#ifdef OPENMP
#include <omp.h>
#endif

#include <string>
#include <cassert>
#include <cmath>
#include "def.h"

namespace sym_lib {
 using namespace std;
 using namespace format;


 // TODO eventually use CSC / CSR class
 // TODO test for these two format change
 int
 csc_to_csr(int nrow, int ncol, int *Ap, int *Ai, double *Ax, int *&rowptr,
            int *&colind, double *&values) {
  // count row entries to generate row ptr
  int nnz = Ap[ncol];
  int *rowCnt = new int[nrow]();
  for (int i = 0; i < nnz; i++)
   rowCnt[Ai[i]]++;

  rowptr = new int[nrow + 1]();
  int counter = 0;
  for (int i = 0; i < nrow; i++) {
   rowptr[i] = counter;
   counter += rowCnt[i];
  }
  rowptr[nrow] = nnz;

  colind = new int[nnz]();
  values = new double[nnz]();

  std::fill_n(rowCnt, nrow, 0);
  //memset(rowCnt, 0, sizeof(int) * nrow);
  for (int i = 0; i < ncol; i++) {
   for (int j = Ap[i]; j < Ap[i + 1]; j++) {
    int row = Ai[j];
    int index = rowptr[row] + rowCnt[row];
    colind[index] = i;
    values[index] = Ax[j];
    rowCnt[row]++;
   }
  }
  delete[]rowCnt;

  return 0;
 }

 CSR* csc_to_csr(CSC* A) {
  // count row entries to generate row ptr
  int nnz = A->p[A->n];
  int *rowCnt = new int[A->n]();
  for (int i = 0; i < nnz; i++)
   rowCnt[A->i[i]]++;

  CSR *B = new CSR(A->n,A->m,A->nnz,A->is_pattern);
  int *rowptr = B->p; //new int[nrow + 1]();
  size_t ncol = B->n;
  size_t nrow = B->m;
  int counter = 0;
  for (int i = 0; i < (int)nrow; i++) {
   rowptr[i] = counter;
   counter += rowCnt[i];
  }
  rowptr[nrow] = nnz;

  int *colind = B->i;
  double *values = B->x;
  std::fill_n(rowCnt, nrow, 0);
  //memset(rowCnt, 0, sizeof(int) * nrow);
  for (int i = 0; i < (int)ncol; i++) {
   for (int j = A->p[i]; j < A->p[i + 1]; j++) {
    int row = A->i[j];
    int index = rowptr[row] + rowCnt[row];
    colind[index] = i;
    if(!B->is_pattern)
     values[index] = A->x[j];
    rowCnt[row]++;
   }
  }
  delete[]rowCnt;
  return B;
 }

 int csr_to_csc(int nrow, int ncol, int *Ai, int *Ap, double *Ax, int *&colptr,
                int *&rowind, double *&values) {
  return csc_to_csr(ncol, nrow, Ap, Ai, Ax, colptr, rowind, values);
 }

 CSC* csr_to_csc(CSR *A) {
  CSC *B = new CSC(A->m,A->n,A->nnz,NULLPNTR,NULLPNTR,A->stype);
  B->pre_alloc = false;
  B->is_pattern = false;
  int st = csc_to_csr(A->m, A->n, A->p, A->i, A->x, B->p, B->i, B->x);
  B->nnz = B->p[B->n];
  return B;
 }

 CSC *tridiag(int n, double a0, double a1, double a2) {
  CSC *A = new CSC(n, n, 3 * n - 2);
  A->is_pattern = true;
  A->stype = -1;

  int ind = 0;
  A->p[0] = 0;
  for(int i = 0; i < n; i++) {
   if(i == 0 || i == n-1) {
    A->p[i+1] = A->p[i] + 2;
    if(i == 0) {
     A->i[ind] = i;
     A->x[ind] = a1;
     A->i[ind+1] = i+1;
     A->x[ind+1] = a2;
    } else {
     A->i[ind] = i-1;
     A->x[ind] = a0;
     A->i[ind+1] = i;
     A->x[ind+1] = a1;
    }
    ind += 2;
   } else {
    A->p[i+1] = A->p[i] + 3;
    A->i[ind + 0] = i-1;
    A->i[ind + 1] = i;
    A->i[ind + 2] = i+1;
    A->x[ind + 0] = a0;
    A->x[ind + 1] = a1;
    A->x[ind + 2] = a2;
    ind += 3;
   }
  }

  return A;
 }


/* CSC* tree_to_csc(int n, int *tree){
  int *Ai = new int[n](), *Ap = new int[n+1](), *cnt = new int[n]();
  populate_children(n, tree, Ap, Ai, cnt);
  CSC *T = new CSC(n,n,n-1,Ap,Ai,0);
  delete []cnt;
  return T;
 }*/

 CSC* transpose_general(CSC *A){
  if(!A)
   return NULLPNTR;
  size_t row = A->m;
  size_t col = A->n;
  size_t  colT = row;
  CSC *AT=NULLPNTR;
  if(row == 0 || col==0){
   return AT;
  }
  AT = new CSC(col,row,A->nnz,A->is_pattern);
  int *col_cnt = new int[row]();//cols will be rows
  for (int i = 0; i < A->p[col]; ++i) {
   col_cnt[A->i[i]]++;
  }
  AT->p[0]=0;
  //Determining column pointer of AT
  for (int j = 1; j < (int)colT+1; ++j) {
   AT->p[j] = AT->p[j-1] + col_cnt[j-1];
  }
  std::fill_n(col_cnt,colT,0);
  //Determining row pointer of AT
  for (int k = 0; k < (int)col; ++k) {
   for (int i = A->p[k]; i < A->p[k+1]; ++i) {
    int beg = A->i[i];
    assert(AT->p[beg]+col_cnt[beg] < A->nnz);
    AT->i[AT->p[beg]+col_cnt[beg]] = k;
    if(!AT->is_pattern)
     AT->x[AT->p[beg]+col_cnt[beg]] = A->x[i];
    col_cnt[beg]++;
   }
  }
  AT->stype = A->stype;
  delete []col_cnt;
  return AT;
 }


 CSC* make_half(size_t An, int *Ap, int *Ai, double *Ax, bool lower){
  size_t Bn = An;
  size_t Bnz = 0;
  for (int i = 0; i < (int)An; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    bool l_u = lower ? Ai[j]>=i : Ai[j]<=i;
    if(l_u){
     Bnz ++;
    }
   }
  }
  CSC *B = new CSC(Bn,Bn,Bnz,Ax==NULLPNTR);
  int *Bp = B->p;
  int *Bi = B->i;
  double *Bx = B->x;
  int nnz_cnt=0;
  for (int i = 0; i < An; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    bool l_u = lower ? Ai[j]>=i : Ai[j]<=i;
    if(l_u){
     Bi[nnz_cnt] = Ai[j];
     if(Ax!=NULLPNTR) Bx[nnz_cnt] = Ax[j];
     nnz_cnt++;
    }
   }
   Bp[i+1] = nnz_cnt;
  }
  B->stype = -1;
  return B;
 }

 CSC *make_full(CSC *A) {
  if(A->stype == format::GENERAL) {
   std::cerr << "Not symmetric\n";
   return nullptr;
  }

  CSC *Afull = new CSC(A->m, A->n, A->nnz * 2 - A->n);
  auto ind = new int[A->n]();

  for(size_t i = 0; i < A->n; i++) {
   for(size_t p = A->p[i]; p < A->p[i+1]; p++) {
    int row = A->i[p];
    ind[i]++;
    if(row != i)
     ind[row]++;
   }
  }
  Afull->p[0] = 0;
  for(size_t i = 0; i < A->n; i++)
   Afull->p[i+1] = Afull->p[i] + ind[i];

  for(size_t i = 0; i < A->n; i++)
   ind[i] = 0;
  for(size_t i = 0; i < A->n; i++) {
   for(size_t p = A->p[i]; p < A->p[i+1]; p++) {
    int row = A->i[p];
    int index = Afull->p[i] + ind[i];
    Afull->i[index] = row;
    Afull->x[index] = A->x[p];
    ind[i]++;
    if(row != i) {
     index = Afull->p[row] + ind[row];
     Afull->i[index] = i;
     Afull->x[index] = A->x[p];
     ind[row]++;
    }
   }
  }
  delete[]ind;

  return Afull;
 }


 CSC *transpose_symmetric(CSC *A, int *Perm) {
  int *Ap, *Ai, *Fi, *Fp;
  double *Ax, *Fx;
  int  *Wi, *Pinv, *Iwork;
  int p, pend, upper, permute, jold, i, j, k, iold, fp;
  size_t n;
  size_t s;
  CSC *F;
  int stype;
  stype = A->stype;
  int is_np = !A->is_pattern;
  if (stype == 0 || A->m != A->n) {
   return NULLPNTR;
  }
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  n = A->m;
  permute = (Perm != NULLPNTR);
  upper = (A->stype > 0);
  F = new CSC(A->n,A->m,A->nnz);
  Fp = F->p;        /* size A->nrow+1, row pointers of F */
  Fi = F->i;
  Fx = F->x;
  s = n +  ((Perm != NULLPNTR) ? n : 0) ;
  Iwork = new int[s]();
  Wi = Iwork;
  Pinv = Iwork + n; // unused if Perm NULL

  //check Perm and construct inverse permutation
  if (permute) {
   for (i = 0; i < n; i++) {
    Pinv[i] = EMPTY;
   }
   for (k = 0; k < n; k++) {
    i = Perm[k];
    if (i < 0 || i > n || Pinv[i] != EMPTY) {
     return NULLPNTR;
    }
    Pinv[i] = k;
   }
  }

  // count the entries in each row of F
  for (i = 0; i < n; i++) {
   Wi[i] = 0;
  }
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
       Wi[std::min(i, j)]++;
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
       Wi[std::max(i, j)]++;
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
  // compute the row pointers
  p = 0;
  for (i = 0; i < n; i++) {
   Fp[i] = p;
   p += Wi[i];
  }
  Fp[n] = p;
  for (i = 0; i < n; i++) {
   Wi[i] = Fp[i];
  }
  if (p > F->nnz) {
   return NULLPNTR;
  }

/// Transpose the matrix
  if (permute) {
   if (upper) {
    /* permuted, upper */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     p = Ap[jold];
     pend = Ap[jold + 1];
     for (; p < pend; p++) {
      iold = Ai[p];
      if (iold <= jold) {
       i = Pinv[iold];
       if (i < j) {
        fp = Wi[i]++;
        Fi[fp] = j;
        if(is_np)
         Fx[fp] = Ax[p];
       } else {
        fp = Wi[j]++;
        Fi[fp] = i;
        if(is_np)
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
     pend = Ap[jold + 1];
     for (; p < pend; p++) {
      iold = Ai[p];
      if (iold >= jold) {
       i = Pinv[iold];
       if (i > j) {
        fp = Wi[i]++;
        Fi[fp] = j;
        if(is_np)
         Fx[fp] = Ax[p];
       } else {
        fp = Wi[j]++;
        Fi[fp] = i;
        if(is_np)
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
     pend = Ap[j + 1];
     for (; p < pend; p++) {
      i = Ai[p];
      if (i <= j) {
       fp = Wi[i]++;
       Fi[fp] = j;
       if(is_np)
        Fx[fp] = Ax[p];
      }
     }
    }
   } else {
    /* unpermuted, lower */
    for (j = 0; j < n; j++) {
     p = Ap[j];
     pend = Ap[j + 1];
     for (; p < pend; p++) {
      i = Ai[p];
      if (i >= j) {
       fp = Wi[i]++;
       Fi[fp] = j;
       if(is_np)
        Fx[fp] = Ax[p];
      }
     }
    }
   }
  }
  F->stype = -A->stype; // flip the stype
  delete []Iwork;
  return F;
 }

 double residual(int beg_idx, int end_idx, const double *x0, const double *x1) {
  double max = 0.0;
  for (int i = beg_idx; i < end_idx; i++)
   if (std::fabs(x0[i] - x1[i]) > max)
    max = std::fabs(x0[i] - x1[i]);
  return max;
 }

 // FIXME: different types of norms
 double norm(int size, double *vec) {
  double max = 0.0;
  for(int i = 0; i < size; i++) {
   if(std::fabs(vec[i]) > max)
    max = std::fabs(vec[i]);
  }
  return max;
 }

 CSC *copy_sparse(CSC *A){
  if(!A)
   return NULLPNTR;
  if(A->n ==0 || A->m ==0)
   return NULLPNTR;
  CSC *clone = new CSC(A->m,A->n,A->nnz);
  for (int i = 0; i < A->n+1; ++i) {
   clone->p[i] = A->p[i];
  }
  for (int j = 0; j < A->nnz; ++j) {
   clone->i[j] = A->i[j];
   clone->x[j] = A->x[j];
  }
  clone->stype = A->stype;
  return clone;
 }

 Dense *copy_dense(Dense *A){
  if(!A)
   return NULLPNTR;
  if(A->row ==0 || A->col == 0)
   return NULLPNTR;
  auto *clone = new Dense(A->row, A->col, A->lda);
  for (int i = 0; i < A->col * A->row; ++i) {
   clone->a[i] = A->a[i];
  }
  return clone;
 }

 void copy_from_to(CSR *src, CSR *dst){
  for (int i = 0; i < src->n; ++i) {
   dst->p[i] = dst->p[i];
   for (int j = src->p[i]; j < src->p[i + 1]; ++j) {
    dst->i[j] = src->i[j];
    dst->x[j] = src->x[j];
   }
  }
 }


 void copy_vector_dense(size_t beg, size_t end, const double *src, double *dst) {
  for(size_t i = beg; i < end; i++)
   dst[i] = src[i];
 }

 template<class type>  type *copy_array(int s, type *val){
  type *clone = new type[s];
  for (int i = 0; i < s; ++i) {
   clone[i] = val[i];
  }
  return clone;
 }

 template<class type>  type *concat_arrays(int s1, type *val1, int s2, type *val2){
  type *clone = new type[s1+s2];
  for (int i = 0; i < s1; ++i) {
   clone[i] = val1[i];
  }
  for (int j = s1; j < s1 + s2; ++j) {
   clone[j] = val2[j];
  }
  return clone;
 }

 int number_empty_col(CSC *A){
  int cnt=0;
  for (int i = 0; i < A->n; ++i) {
   if(A->p[i+1]-A->p[i] == 0){
    cnt++;
   }
  }
  return cnt;
 }




 void sparse2dense(CSC *A, double *D){
  std::fill_n(D,A->n*A->m,0);
  for (int i = 0; i < A->n; ++i) {
   for (int j = A->p[i]; j < A->p[i+1]; ++j) {
    int r = A->i[j];
    D[i*A->m+r] = A->x[j];
    if(A->stype==-1 || A->stype==1){
     D[r*A->n+i] = A->x[j];
    }
   }
  }
 }


 CSC *diagonal(int n, double val){
  CSC *A = new CSC(n,n,n);
  A->p[0] = 0;
  for (int i = 0; i < n; ++i) {
   A->i[i] = i;
   A->p[i+1] = A->p[i] + 1;
   A->x[i] = val;
  }
  assert(A->p[n] == n);
  return A;
 }



 void compute_nnz_per_col(CSC *A, double *nnz_cnt){
  for (int i = 0; i < A->n; ++i) {
   nnz_cnt[i] = A->p[i+1] - A->p[i];
  }
 }

/*
 * Concatenates to non-square matrices in CSC and generates
 * a new non-square CSC matrix. A|B = [A]
 *                                    [B]
 */
 CSC* concatenate_two_CSC(CSC *A, CSC *B){
  CSC *SM = NULLPNTR;
  if(!A && !B){
   return NULLPNTR;
  }else if(!A){
   SM = copy_sparse(B);
   return SM;
  } else if(!B){
   SM = copy_sparse(A);
   return SM;
  }
  size_t SM_nz;
  int *SMp;
  int *SMi;
  double *SMx;
  size_t SM_size=0;
  if(B->m>0){
   SM_size = A->n; //cols
   SMp = new int[SM_size+1];
   SMp[0] = 0;
   for (int i = 1; i < SM_size+1; ++i) {
    SMp[i] = SMp[i-1] + (A->p[i] - A->p[i-1]) +
             (B->p[i] - B->p[i-1]);
   }
  }
  SM_nz = SMp[SM_size];
  SMi = new int[SM_nz];
  SMx = new double[SM_nz]();

  int base1=A->m;
  size_t stp = 0;
  for (int j = 0; j < SM_size; ++j) {
   stp = SMp[j];
   //Adding A
   for (int i = A->p[j]; i < A->p[j+1]; ++i) {
    SMi[stp] = A->i[i];
    if(!A->is_pattern)
     SMx[stp] = A->x[i];
    stp++;
   }
   base1=A->m;
   //Adding B values
   if(B->m > 0){
    for (int i = B->p[j]; i < B->p[j+1]; ++i) {
     SMi[stp] = base1 + B->i[i];
     if(!B->is_pattern)
      SMx[stp] = B->x[i];
     stp++;
    }
   }
   assert(stp == SMp[j+1]);
  }
  SM = new CSC(A->m+B->m, SM_size, SM_nz, SMp, SMi, SMx);
  SM->is_pattern= false;
  SM->stype = A->stype;
  return SM;
 }

 Dense* concatenate_two_dense(Dense *a, Dense *b){
  Dense *ret;
  if(!a && !b)
   return NULLPNTR;
  if(!a){
   ret = copy_dense(b);
   return ret;
  }
  if(!b){
   ret = copy_dense(a);
   return ret;
  }
  assert(a->col = b->col);
  int nrows = a->row + b->row;
  ret = new Dense(nrows, a->col, 1);
  for (int i = 0; i < a->col; ++i) {
   for (int j = 0; j < a->row; ++j) {
    ret->a[i*nrows + j] = a->a[j];
   }
   for (int k = 0; k < b->row; ++k) {
    ret->a[i*nrows + a->row +k] = b->a[k];
   }
  }
  return ret;
 }

 template <class T> bool are_equal(T *m1, T *m2){
  if( m1 == NULLPNTR && m2 == NULLPNTR)
   return true;
  return m1 ? m1->equality_check(m2) : m2->equality_check(m1);
 }

}

