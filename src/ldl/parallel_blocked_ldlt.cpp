//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/parallel_blocked_ldlt.h"

#include <cassert>
#include <chrono>
#include <nasoq/common/Sym_BLAS.h>




#include "nasoq/common/Reach.h"

namespace nasoq {

 bool
 ldl_left_sn_parallel_01(int n, int *c, int *r, double *values, size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                         double *D, int *blockSet, int supNo, double *timing, int *aTree, int *cT, int *rT,
                         int *col2Sup, int nLevels, int *levelPtr, int *levelSet, int nPar, int *parPtr, int *partition,
                         int chunk, int threads, int super_max, int col_max, int &nbpivot, double threshold) {
  /*
   * For timing using BLAS
   */
  const int incx = 1;
  int top = 0;
  int *xi = new int[2 * supNo]();
  //int super_max = 64; //tunig parameter for the max size of supernodes TODO: find it in analysis
  //int col_max = n;
  int *map; //= new int[n]();
  double *contribs; //= new double[super_max*col_max]();
  double *trn_diag; //= new double[super_max*col_max]();
  int info;
  double one[2], zero[2];
  one[0] = 1.0;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
  one[1] = 0.;
  zero[0] = 0.;     /* BETA for *syrk, *herk, and *gemm */
  zero[1] = 0.;
  int *ipiv = new int[n]();
  std::chrono::time_point<std::chrono::system_clock> start, end, startin, endin;
  std::chrono::duration<double> elapsed_seconds;
  double duration4 = 0, duration3 = 0, duration2 = 0, duration1 = 0;
#ifdef TIMING
  start = std::chrono::system_clock::now();
#endif
  for (int i1 = 0; i1 < nLevels - 1; ++i1) {
#pragma omp parallel //shared(lValues)//private(map, contribs)
   {
#pragma omp  for schedule(dynamic) private(map, trn_diag, contribs, xi, startin, endin, duration2)
    for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
#ifdef BLASTIMING
     int threadID = omp_get_thread_num();
     std::chrono::time_point<std::chrono::system_clock> startBlas, endBlas;
#endif
     map = new int[n]();
     contribs = new double[super_max * col_max]();
     xi = new int[2 * supNo]();
     trn_diag = new double[super_max * col_max]();
     //int pls = levelSet[j1];
#ifdef TIMING1
     startin = std::chrono::system_clock::now();
#endif
//#pragma omp parallel for schedule(static,chunk)private(thth)
     for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
      int s = partition[k1] + 1;

      int curCol = s != 0 ? blockSet[s - 1] : 0;
      int nxtCol = blockSet[s];
      int supWdt = nxtCol - curCol;
      int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
      for (int i = Li_ptr[curCol], cnt = 0; i < Li_ptr[nxtCol]; ++i) {
       map[lR[i]] = cnt++;//mapping L rows position to actual row idx
      }

      //copy the columns from A to L
      for (int i = curCol; i < nxtCol; ++i) {//Copy A to L
       int pad = i - curCol;
       for (int j = c[i]; j < c[i + 1]; ++j) {
        // if(r[j]>=i)//does not need to save upper part.
        lValues[lC[i] + map[r[j]]] = values[j];
        //   else
        //      printf("dddd\n");
       }
      }
#if DEBUG
      top = ereach_sn(supNo,c,r,curCol,nxtCol,col2sup, eTree,xi,xi+supNo);
            if(supNo-top != prunePtr[s]-prunePtr[s-1])
                printf("sss");
#endif
      double *src, *cur = &lValues[lC[curCol]];//pointing to first element of the current supernode

//#ifndef PRUNE
      top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
      assert(top >= 0);
      //if(s==2){top =2; xi[top] = 0;}
      for (int i = top; i < supNo; ++i) {
       int lSN = xi[i];
/*#else
      for (int i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
         int lSN = pruneSet[i];
#endif*/
       int nSupRs = 0;
#if DEBUG
       if(xi[top++] != lSN)
                    printf("fail");
#endif
       int cSN = blockSet[lSN];//first col of current SN
       int cNSN = blockSet[lSN + 1];//first col of Next SN
       int Li_ptr_cNSN = Li_ptr[cNSN];
       int Li_ptr_cSN = Li_ptr[cSN];
       int nSNRCur = Li_ptr_cNSN - Li_ptr_cSN;
       int supWdts = cNSN - cSN;//The width of current src SN
       int lb = 0, ub = 0;
       bool sw = true;
       int beg_col = cSN, end_col = 0;
       for (int j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
        //finding the overlap between curCol and curCol+supWdt in the src col
        if (lR[j] >= curCol && sw) {
         //src*transpose(row lR[j])
         lb = j - Li_ptr_cSN;
         sw = false;
        }
        if (lR[j] < curCol + supWdt && !sw) {
         ub = j - Li_ptr_cSN;
        }
       }
       nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
       int ndrow1 = ub - lb + 1;
       int ndrow3 = nSupRs - ndrow1;
       src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
       double *srcL = &lValues[lC[cSN] + ub + 1];
       // multiplying L * D
       for (int l = 0; l < supWdts; ++l) {
        double tmp = D[cSN + l];
        for (int l1 = 0; l1 < nSupRs; ++l1) {
         trn_diag[l * nSupRs + l1] = tmp * src[l * nSNRCur + l1];
        }
       }
       /*dgemm("N", "C", &ndrow3, &ndrow1, &supWdts, one, srcL, &nSNRCur,
             src, &nSNRCur, zero, &contribs[ndrow1], &nSupRs);*/
#ifdef OPENBLAS
       cblas_dgemm(CblasColMajor,CblasNoTrans,CblasConjTrans, nSupRs, ndrow1, supWdts, 1.0, trn_diag, nSupRs,
             src, nSNRCur, 0.0, contribs, nSupRs);
#else
       SYM_DGEMM("N", "C", &nSupRs, &ndrow1, &supWdts, one, trn_diag, &nSupRs,
             src, &nSNRCur, zero, contribs, &nSupRs);
#endif

       //copying contrib to L
       for (int i = 0; i < ndrow1; ++i) {//Copy contribs to L
        int col = map[lR[Li_ptr_cSN + i + lb]];//col in the SN
        //double ddiag = 1.0 ;/// D[col];
        for (int j = i; j < nSupRs; ++j) {
         int cRow = lR[Li_ptr_cSN + j + lb];//corresponding row in SN
         //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
         cur[col * nSupR + map[cRow]] -= contribs[i * nSupRs + j];
        }
       }
      }
      sym_sytrf(cur, supWdt, nSupR, &nbpivot, threshold);

      // Making L*D
      int rowNo = nSupR - supWdt;
      for (int l = 0; l < supWdt; ++l) {
       double tmp = cur[l + l * nSupR];
       D[curCol + l] = tmp;
       double *stCol = trn_diag + l * supWdt + l;
       double *curCol = cur + l * nSupR + l;
       *stCol = tmp;
       for (int l1 = 0; l1 < supWdt - l - 1; ++l1) {
        *(++stCol) = tmp * *(++curCol);
       }
      }
#ifdef OPENBLAS
      cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans, CblasNonUnit, rowNo, supWdt, 1.0,
                trn_diag, supWdt, &cur[supWdt], nSupR);
#else
      SYM_DTRSM("R", "L", "C", "N", &rowNo, &supWdt, one,
            trn_diag, &supWdt, &cur[supWdt], &nSupR);
#endif

      for (int k = 0; k < supWdt; ++k) {
       cur[k * nSupR + k] = 1.0;
      }
      //copying 1/Di into D
      /*for (int l = 0; l < supWdt; ++l) {
       D[curCol+l] = one[0] / cur[l + l*nSupR];
      }*/


/*  for (int k = 0; k < nSupR * supWdt; ++k) {
   std::cout<<cur[k]<<",";
  }
  std::cout<<"==== \n";*/
     }
     delete[]contribs;
     delete[]trn_diag;
     delete[]xi;
     delete[]map;
    }
#ifdef TIMING1
    endin = std::chrono::system_clock::now();
    elapsed_seconds = endin-startin;
    duration1=elapsed_seconds.count();
    int thth2=omp_get_thread_num();
    std::cout<<"**"<<thth2<<" : "<<j1<<" "<<duration1<<"\n";
#endif

   }
  }

#if 1
  //LAst iteration
#ifdef OPENBLAS
  openblas_set_num_threads(threads);
#else
  MKL_Domain_Set_Num_Threads(threads, MKL_DOMAIN_BLAS);
#endif


  map = new int[n]();
  contribs = new double[super_max * col_max]();
  xi = new int[2 * supNo]();
  trn_diag = new double[super_max * col_max]();
  for (int j1 = levelPtr[nLevels - 1]; j1 < levelPtr[nLevels]; ++j1) {
#ifdef TLAST
   start = std::chrono::system_clock::now();
#endif
   for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
    int s = partition[k1] + 1;

    int curCol = s != 0 ? blockSet[s - 1] : 0;
    int nxtCol = blockSet[s];
    int supWdt = nxtCol - curCol;
    int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
    for (int i = Li_ptr[curCol], cnt = 0; i < Li_ptr[nxtCol]; ++i) {
     map[lR[i]] = cnt++;//mapping L rows position to actual row idx
    }

    //copy the columns from A to L
    for (int i = curCol; i < nxtCol; ++i) {//Copy A to L
     int pad = i - curCol;
     for (int j = c[i]; j < c[i + 1]; ++j) {
      // if(r[j]>=i)//does not need to save upper part.
      lValues[lC[i] + map[r[j]]] = values[j];
      //   else
      //      printf("dddd\n");
     }
    }
#if DEBUG
    top = ereach_sn(supNo,c,r,curCol,nxtCol,col2sup, eTree,xi,xi+supNo);
            if(supNo-top != prunePtr[s]-prunePtr[s-1])
                printf("sss");
#endif
    double *src, *cur = &lValues[lC[curCol]];//pointing to first element of the current supernode

//#ifndef PRUNE
    top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
    assert(top >= 0);
    //if(s==2){top =2; xi[top] = 0;}
    for (int i = top; i < supNo; ++i) {
     int lSN = xi[i];
/*#else
      for (int i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
         int lSN = pruneSet[i];
#endif*/
     int nSupRs = 0;
#if DEBUG
     if(xi[top++] != lSN)
                    printf("fail");
#endif
     int cSN = blockSet[lSN];//first col of current SN
     int cNSN = blockSet[lSN + 1];//first col of Next SN
     int Li_ptr_cNSN = Li_ptr[cNSN];
     int Li_ptr_cSN = Li_ptr[cSN];
     int nSNRCur = Li_ptr_cNSN - Li_ptr_cSN;
     int supWdts = cNSN - cSN;//The width of current src SN
     int lb = 0, ub = 0;
     bool sw = true;
     int beg_col = cSN, end_col = 0;
     for (int j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
      //finding the overlap between curCol and curCol+supWdt in the src col
      if (lR[j] >= curCol && sw) {
       //src*transpose(row lR[j])
       lb = j - Li_ptr_cSN;
       sw = false;
      }
      if (lR[j] < curCol + supWdt && !sw) {
       ub = j - Li_ptr_cSN;
      }
     }
     nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
     int ndrow1 = ub - lb + 1;
     int ndrow3 = nSupRs - ndrow1;
     src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
     double *srcL = &lValues[lC[cSN] + ub + 1];
     // multiplying L * D
     for (int l = 0; l < supWdts; ++l) {
      double tmp = D[cSN + l];
/*     double *dst_tmp = &rn_diag[l * nSupRs];
     double *src_tmp = &src[l * nSNRCur];
     dscal(nSupRs, )*/
      for (int l1 = 0; l1 < nSupRs; ++l1) {
       trn_diag[l * nSupRs + l1] = tmp * src[l * nSNRCur + l1];
      }
     }
     /*dgemm("N", "C", &ndrow3, &ndrow1, &supWdts, one, srcL, &nSNRCur,
           src, &nSNRCur, zero, &contribs[ndrow1], &nSupRs);*/
#ifdef OPENBLAS
     cblas_dgemm(CblasColMajor,CblasNoTrans,CblasConjTrans, nSupRs, ndrow1, supWdts, 1.0, trn_diag, nSupRs,
                 src, nSNRCur, 0.0, contribs, nSupRs);
#else
     SYM_DGEMM("N", "C", &nSupRs, &ndrow1, &supWdts, one, trn_diag, &nSupRs,
           src, &nSNRCur, zero, contribs, &nSupRs);

#endif


     //copying contrib to L
     for (int i = 0; i < ndrow1; ++i) {//Copy contribs to L
      int col = map[lR[Li_ptr_cSN + i + lb]];//col in the SN
      //double ddiag = 1.0 ;/// D[col];
      for (int j = i; j < nSupRs; ++j) {
       int cRow = lR[Li_ptr_cSN + j + lb];//corresponding row in SN
       //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
       cur[col * nSupR + map[cRow]] -= contribs[i * nSupRs + j];
      }
     }
    }
    sym_sytrf(cur, supWdt, nSupR, &nbpivot, threshold);

    // Making L*D
    int rowNo = nSupR - supWdt;
    for (int l = 0; l < supWdt; ++l) {
     double tmp = cur[l + l * nSupR];
     D[curCol + l] = tmp;
     double *stCol = trn_diag + l * supWdt + l;
     double *curCol = cur + l * nSupR + l;
     *stCol = tmp;
     for (int l1 = 0; l1 < supWdt - l - 1; ++l1) {
      *(++stCol) = tmp * *(++curCol);
     }
    }
#ifdef OPENBLAS
    cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans, CblasNonUnit, rowNo, supWdt, 1.0,
                trn_diag, supWdt, &cur[supWdt], nSupR);
#else
    SYM_DTRSM("R", "L", "C", "N", &rowNo, &supWdt, one,
          trn_diag, &supWdt, &cur[supWdt], &nSupR);

#endif

    for (int k = 0; k < supWdt; ++k) {
     cur[k * nSupR + k] = 1.0;
    }
    //copying 1/Di into D
    /*for (int l = 0; l < supWdt; ++l) {
     D[curCol+l] = one[0] / cur[l + l*nSupR];
    }*/


/*  for (int k = 0; k < nSupR * supWdt; ++k) {
   std::cout<<cur[k]<<",";
  }
  std::cout<<"==== \n";*/
   }
#ifdef TLAST
   end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  duration1=elapsed_seconds.count();
  std::cout<<"++ " <<duration1<<"\n";
#endif
  }
#endif
  delete[]contribs;
  delete[]trn_diag;
  delete[]xi;
  delete[]map;
  return true;
 }
}