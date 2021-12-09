//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/ldlt_check.h"

#include "nasoq/common/def.h"
#include "nasoq/common/Transpose.h"
#include "nasoq/matrixMatrix/spmm.h"

namespace nasoq {

 bool ldlt_check(int n, size_t *Lp, int *Li, double *valL, double *D, size_t *Ap, int *Ai, double *valA) {
  double *ll_col = new double[n]();
  double *a_col = new double[n]();
  double *row_comp = new double[n]();
  double *valTmp = new double[Lp[n]];
  int status = 0;

  CSC *L_csc = new CSC;
  L_csc->nzmax = Lp[n];
  L_csc->ncol = L_csc->nrow = n;
  L_csc->stype = -1;
  L_csc->xtype = CHOLMOD_REAL;
  L_csc->packed = TRUE;
  L_csc->p = new int[n + 1];
  for (int m = 0; m <= n; ++m) {//FIXME
   L_csc->p[m] = static_cast<int >(Lp[m]);
  }
  L_csc->i = Li;
  L_csc->x = valL;
  L_csc->nz = NULL;
  L_csc->sorted = TRUE;
  CSC *LT = ptranspose(L_csc, 2, NULL, NULL, 0, status);
  // sparse matrix times diagonal
  spmdm(n, Lp, Li, valL, D, valTmp);
  for (int i = 0; i < n; ++i) {

   // col i of L'
   for (int j = LT->p[i]; j < LT->p[i + 1]; ++j) {
    ll_col[LT->i[j]] = LT->x[j];
   }
   // col i of A
   for (int k = Ap[i]; k < Ap[i + 1]; ++k) {
    a_col[Ai[k]] = valA[k];
   }
   // LD time L'[i,:] = col i of LDL'
   spmv_csc(n, Lp, Li, valTmp, ll_col, row_comp);

   //Verifying
   for (int l = i; l < n; ++l) {
    // Warning: Sometimes the matrix condition number causes lose of accuracy
    // not necessarily the implementation.
    if (row_comp[l] - a_col[l] > 1e-3) {
     printf("row %d : col %d \n", i, l);
     return false;
    }
   }
   //reseting tmp arrays.
   for (int j = LT->p[i]; j < LT->p[i + 1]; ++j) {
    ll_col[LT->i[j]] = 0.0;
   }
   // col i of A
   for (int k = Ap[i]; k < Ap[i + 1]; ++k) {
    a_col[Ai[k]] = 0.0;
   }
   std::fill(row_comp, row_comp + n, 0);

  }
  delete[]ll_col;
  delete[]a_col;
  delete[]row_comp;
  delete[]valTmp;
  delete[]L_csc->p;
  delete L_csc;
  allocateAC(LT, 0, 0, 0, FALSE);
  return true;

 }
}