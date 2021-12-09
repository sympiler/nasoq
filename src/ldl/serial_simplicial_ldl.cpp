//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/serial_simplicial_ldl.h"


#include <algorithm>
#include <cassert>

#include "nasoq/common/Reach.h"

namespace nasoq {

 int ldl_left_simplicial_01(int n, int *c, int *r, double *values, int *cT, int *rT, int *lC, int *lR, double *&lValues,
                            double *d, int *eTree, double *ws, int *ws_int) {
  //ws n, ws_int size of 3*n
  /*
   * Performs a Cholesky decomposition on a given matrix (c, r, values), i.e.
   * stored in compressed column format, and produces L, which are
   * stored in column compressed format.
   * (n, c, r, values) : IN : input matrix
   * (lC, lR) : IN : The column and rwo sparsity patterns of L
   * (lValues) : OUT : the nonzero values of the L factor
   * (pruneSet, prunePtr) : IN : the row sparsity pattern of the factor L
   */
  int top = 0;
  double *f = ws;
  int *finger = ws_int;
  int *xi = ws_int + n;
  std::fill_n(f, n, 0);
  std::fill_n(xi, 2 * n, 0);
  //Determining the column pointer
  int spCol = 0;
  for (int colNo = 0; colNo < n; ++colNo) {
   //Uncompress a col into a 1D array
   for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
    f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
   }
#if 0
   for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]-1; ++i) {
   spCol = pruneSet[i];
#endif
   top = ereach(n, cT, rT, colNo, eTree, xi, xi + n);
   //std::cout<<n-top<<";\n";
   for (int i = top; i < n; ++i) {
    spCol = xi[i];
    double facing_val = lValues[finger[spCol]] * d[spCol];
    for (int l = lC[spCol]; l < lC[spCol + 1]; ++l) {
     if (lR[l] > colNo) {
      f[lR[l]] -= lValues[l] * facing_val;
     }
    }
    d[colNo] += facing_val * lValues[finger[spCol]];
    finger[spCol]++;
   }
   finger[colNo]++;//Skip diagonal row
   d[colNo] = f[colNo] - d[colNo];
   double diag = d[colNo];
   //double tmpSqrt = sqrt(f[colNo]);
   f[colNo] = 0;
   lValues[lC[colNo]] = 1;
   for (int j = lC[colNo] + 1; j < lC[colNo + 1]; ++j) {
    lValues[j] = f[lR[j]] / diag;
    f[lR[j]] = 0;
   }
  }
  return 1;
 }

 int ldl_left_simplicial_02(int n, int *c, int *r, double *values, int *cT, int *rT, int *lC, int *lR, double *&lValues,
                            double *d, int *eTree, double *ws, int *ws_int) {
  //ws n, ws_int size of 3*n
  /*
   * Performs a Cholesky decomposition on a given matrix (c, r, values), i.e.
   * stored in compressed column format, and produces L, which are
   * stored in column compressed format.
   * (n, c, r, values) : IN : input matrix
   * (lC, lR) : IN : The column and rwo sparsity patterns of L
   * (lValues) : OUT : the nonzero values of the L factor
   * (pruneSet, prunePtr) : IN : the row sparsity pattern of the factor L
   */
  int top = 0;
  double *f = ws;
  int *finger = ws_int;
  int *xi = ws_int + n;
  std::fill_n(f, n, 0);
  std::fill_n(xi, 2 * n, 0);
  //Determining the column pointer
  for (int k = 0; k < n; ++k) {
   finger[k] = lC[k];
  }
  int spCol = 0;
  for (int colNo = 0; colNo < n; ++colNo) {
   //emptying col colNo in L
   //Uncompress a col into a 1D array
   for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
    f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
   }
#if 0
   for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]-1; ++i) {
    spCol = pruneSet[i];
#endif
   top = ereach(n, cT, rT, colNo, eTree, xi, xi + n);
   //std::cout<<n-top<<";\n";
   for (int i = top; i < n; ++i) {
    spCol = xi[i];
    bool sw = false;
    double facing_val = 0, tmp = 0;
    int facing_idx = -1;
    for (int l = lC[spCol]; l < lC[spCol + 1]; ++l) {
     if (lR[l] == colNo) {
      facing_val = lValues[l];
      tmp = facing_val * d[spCol];
      facing_idx = l;
      break;
     }
    }
    assert(facing_idx >= 0);
    for (int l = facing_idx + 1; l < lC[spCol + 1]; ++l) {
     f[lR[l]] -= lValues[l] * tmp;
    }
    d[colNo] += facing_val * tmp;
   }
   d[colNo] = f[colNo] - d[colNo];
   double diag = d[colNo];
   //double tmpSqrt = sqrt(f[colNo]);
   f[colNo] = 0;
   lValues[lC[colNo]] = 1;
   for (int j = lC[colNo] + 1; j < lC[colNo + 1]; ++j) {
    lValues[j] = f[lR[j]] / diag;
    f[lR[j]] = 0;
   }
/*  for (int j = lC[colNo] + 1; j < lC[colNo+1]; ++j) {
   std::cout<<lValues[j]<<";";
  }
  std::cout<<"\n";*/
  }
  return 1;
 }
}