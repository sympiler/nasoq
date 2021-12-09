//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/symbolic/performanceModel.h"

#include <cassert>

#include "nasoq/common/Reach.h"

namespace nasoq {

 double computeCostColFact(double colNo, double len) {
  double total = 0;
  total = (double) (len * PERF_COPY(colNo) + PERF_PPF(colNo) + PERF_TRSM(colNo, len));
  return (total > 0) ? total : 0;
 }

 double computeCostColFact2(double colNo, double len) {
  double total = 0;
  total = (double) (OPS_PPF(colNo) + OPS_TRSM(colNo, len));
  return (total > 0) ? total : 0;
 }

 double computeCostColFact3(double colNo, double len) {
  double total = 0;
  total = colNo;
  return (total > 0) ? total : 0;
 }

 double computeCostColFact4(double colNo, double len) {
  double total = 0;
  total = colNo * len;
  return (total > 0) ? total : 0;
 }

 double computeContribCost(double dim1, double midDim, double dim3) {
  double total = 0;
  total = (double) (PERF_GEMM(dim3, midDim, dim1));
  return (total > 0) ? total : 0;
 }

 double computeContribCost2(double dim1, double midDim, double dim3) {
  double total = 0;
  total = (double) (OPS_GEMM(dim3, midDim, dim1));
  return (total > 0) ? total : 0;
 }

 double computeContribCost3(double dim1, double midDim, double dim3) {
  double total = midDim;
  return (total > 0) ? total : 0;
 }

 double computeContribCost4(double dim1, double midDim, double dim3) {
  double total = midDim * (dim1 + dim3);
  return (total > 0) ? total : 0;
 }

 double computeCostperCol(int n, int colNo, int *eTree, int *cT, int *rT, int *xi, int &top) {
  //int *xi = new int[2*n]();
  double total = 0.0;
  top = ereach(n, cT, rT, colNo, eTree, xi, xi + n);
  for (int i = top; i < n; ++i) {
   int spCol = xi[i];
   bool sw = false;
   int ub = 0;
   total += cT[spCol + 1] - cT[spCol];
  }
  //std::cout<<total<<",";
  total += cT[colNo + 1] - cT[colNo];
  return total;
 }

 double
 computeCostperBlock(int n, int curCol, int nxtCol, int *eTree, int *cT, int *rT, int *col2Sup, int *blockSet, int *lR,
                     size_t *Li_ptr, int *xi) {
  //int *xi = new int[2*n]();
  int supWdt = nxtCol - curCol;
  size_t Li_ptr_nxt = Li_ptr[nxtCol];
  size_t Li_ptr_cur = Li_ptr[curCol];
  int supLen = Li_ptr_nxt - Li_ptr_cur;
  double total = 0.0;
  int top = ereach_sn(n, cT, rT, curCol, nxtCol, col2Sup, eTree, xi, xi + n);
  for (int i = top; i < n; ++i) {
   int lSN = xi[i];
   int nSupRs = 0;
   int cSN = blockSet[lSN];//first col of current SN
   int cNSN = blockSet[lSN + 1];//first col of Next SN
   size_t Li_ptr_cNSN = Li_ptr[cNSN];
   size_t Li_ptr_cSN = Li_ptr[cSN];
   int supWdts = cNSN - cSN;//The width of current src SN
   int lb = 0, ub = 0;
   bool sw = true;
   for (size_t j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
    //finding the overlap between curCol and curCol+supWdt in the src col
    if (lR[j] >= curCol && sw) {
     //src*transpose(row lR[j])
     lb = j - Li_ptr_cSN;
     sw = false;
    }
    if (lR[j] < curCol + supWdt && !sw) {
     ub = j - Li_ptr_cSN;
    }
    if (lR[j] >= curCol + supWdt)
     break;
   }
   nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
   int ndrow1 = ub - lb + 1;
   int ndrow3 = nSupRs - ndrow1;
   total += computeContribCost4(ndrow1, supWdts, ndrow3);
  }
  //std::cout<<total<<",";
  total += computeCostColFact4(supWdt, supLen);
  return total;
 }

 double computeComCost(int n, int *c, int *r, double *values, size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                       int *blockSet, int supNo, double *timing, int *aTree, int *cT, int *rT, int *col2Sup,
                       int nLevels, int *levelPtr, int *levelSet, int nPar, int *parPtr, int *partition,
                       long &totalIntraContribNNZs, long &totalInterCotribNNZs, int &numberOfIntraCore,
                       int &numberOfInterCore, long &numberOfEdgeCuts, double &avgNNZperThread) {
  int *xi = new int[2 * supNo]();
  int *col2part = new int[supNo]();
  int *col2level = new int[supNo]();
  int numOfThreads = 6;
  double *costPerThread = new double[numOfThreads];
  for (int i = 0; i < numOfThreads; ++i) {
   costPerThread[i] = 0.0;
  }

  for (int k = 0; k < nLevels; ++k) {
   for (int i = levelPtr[k]; i < levelPtr[k + 1]; ++i) {
    for (int j = parPtr[i]; j < parPtr[i + 1]; ++j) {
     int node = partition[j];
     col2part[node] = i;
     col2level[node] = k;
    }
   }
  }

  double total = 0.0, tmp = 0.0;
  for (int col = 1; col <= supNo; ++col) {
   int curCol = col != 0 ? blockSet[col - 1] : 0;
   int nxtCol = blockSet[col];
   int supWdt = nxtCol - curCol;
   size_t Li_ptr_nxt = Li_ptr[nxtCol];
   size_t Li_ptr_cur = Li_ptr[curCol];
   int supLen = Li_ptr_nxt - Li_ptr_cur;
   int colPar = aTree[col];
   int mappedThread = col2part[col] % numOfThreads;
   assert(mappedThread >= 0 && mappedThread < numOfThreads);
   if (colPar >= 0) {
    if (col2level[col] != col2level[colPar]) {
     numberOfEdgeCuts++;
    }
   }
   assert(costPerThread[mappedThread] >= 0);
   int top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
   for (int i = top; i < supNo; ++i) {
    int lSN = xi[i];
    int nSupRs = 0;
    int cSN = blockSet[lSN];//first col of current SN
    int cNSN = blockSet[lSN + 1];//first col of Next SN
    size_t Li_ptr_cNSN = Li_ptr[cNSN];
    size_t Li_ptr_cSN = Li_ptr[cSN];
    int supWdts = cNSN - cSN;//The width of current src SN
    int lb = 0, ub = 0;
    bool sw = true;
    for (size_t j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
     //finding the overlap between curCol and curCol+supWdt in the src col
     if (lR[j] >= curCol && sw) {
      //src*transpose(row lR[j])
      lb = j - Li_ptr_cSN;
      sw = false;
     }
     if (lR[j] < curCol + supWdt && !sw) {
      ub = j - Li_ptr_cSN;
     }
     if (lR[j] >= curCol + supWdt)
      break;
    }
    nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
    int ndrow1 = ub - lb + 1;
    int ndrow3 = nSupRs - ndrow1;
    tmp = computeContribCost4(ndrow1, supWdts, ndrow3);
    total += tmp;
    if ((col2level[curCol] == col2level[lSN] &&
         col2part[curCol] == col2part[lSN])
        || col2part[lSN] == 0) {
     totalIntraContribNNZs += ndrow1 * supWdts;
     numberOfIntraCore++;
    } else {
     totalInterCotribNNZs += ndrow1 * supWdts;
     numberOfInterCore++;
     assert(costPerThread[mappedThread] >= 0);
    }
   }
   //std::cout<<total<<",";
   tmp = computeCostColFact4(supWdt, supLen);
   total += tmp;
   assert(tmp > 0);
   assert(costPerThread[mappedThread] >= 0);
   costPerThread[mappedThread] += tmp;
  }
#if 0 //NOT USED for now
  avgNNZperThread=0;
  int usedThreads =0 ;
  for (int l = 0; l < numOfThreads; ++l) {
   if(costPerThread[l]>0){
    avgNNZperThread+=costPerThread[l];
    usedThreads++;
   }

  }
  avgNNZperThread /= (usedThreads>0 ? usedThreads : 1);
#endif
  delete[]costPerThread;
  delete[]col2level;
  delete[]col2part;
  delete[]xi;
  return total;
 }
}