//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/Parallel_update_simplicial.h"

#include <cassert>

#ifdef OPENMP
#include "omp.h"
#endif


#include "nasoq/common/Reach.h"

namespace nasoq {

 int update_ldl_parallel_left_simplicial_01(int n, int *c, int *r, double *values, int *cT, int *rT, int *lC, int *lR,
                                            double *&lValues, double *d, int *eTree, bool *marked, int nLevels,
                                            int *levelPtr, int nPar, int *parPtr, int *partition) {
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

  double *f;
  int *xi;

#ifdef OPENMP
  omp_set_nested(1);
#endif
  for (int i1 = 0; i1 < nLevels; ++i1) {
#pragma omp parallel private(f, xi)
   {
#pragma omp  for schedule(static) private(f, xi)
    for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
     f = new double[n]();
     xi = new int[2 * n]();
     //int pls = levelSet[j1];
     for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
      int colNo = partition[k1];
      if (!marked[colNo])
       continue;
      int spCol = 0;
      int top = 0;

      //emptying col colNo in L
      for (int k = lC[colNo]; k < lC[colNo + 1]; ++k) {
       lValues[k] = 0;
      }
      d[colNo] = 0;
      //Uncompress a col into a 1D array
      for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
       f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
      }
/*#if 0
     for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]-1; ++i) {
   spCol = pruneSet[i];
#endif*/
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
/*     std::cout<<" -> "<<colNo<<"\n";
       for (int j = lC[colNo] + 1; j < lC[colNo+1]; ++j) {
   std::cout<<lValues[j]<<";";
  }
  std::cout<<"\n";*/
     }
     delete[]f;
     delete[]xi;
    }
   }
  }
  return 1;
 }
}