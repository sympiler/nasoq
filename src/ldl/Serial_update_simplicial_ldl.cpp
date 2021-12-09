//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/Serial_update_simplicial_ldl.h"

#include <algorithm>
#include <cassert>

#include "nasoq/common/Reach.h"

namespace nasoq {

 int update_ldl_left_simplicial_01(int n, int *c, int *r, double *values, int *cT, int *rT, int *lC, int *lR,
                                   double *&lValues, double *d, int *eTree, std::vector<int> mod_indices, double *ws,
                                   int *ws_int) {
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
  int colNo = 0;
  //std::cout<<"total size: "<<mod_indices.size();
  for (int cm = 0; cm < mod_indices.size(); ++cm) {
   colNo = mod_indices[cm];
   //emptying col colNo in L
   for (int k = lC[colNo]; k < lC[colNo + 1]; ++k) {
    lValues[k] = 0;
   }
   d[colNo] = 0;
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
/*  std::cout<<" -> "<<colNo<<"\n";
  for (int j = lC[colNo] + 1; j < lC[colNo+1]; ++j) {
   std::cout<<lValues[j]<<";";
  }
  std::cout<<"\n";*/
  }
  return 1;
 }
}