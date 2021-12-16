//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/Serial_update_ldl_static.h"

#include <cassert>
#include <nasoq/common/Sym_BLAS.h>




#include "nasoq/common/Reach.h"

namespace nasoq {

 bool update_ldl_left_sn_01(int n, int *c, int *r, double *values, size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                            double *D, int *blockSet, int supNo, double *timing, int *aTree, int *cT, int *rT,
                            int *col2Sup, std::vector<int> mod_indices, int super_max, int col_max, int &nbpivot,
                            double threshold) {
  /*
   * For timing using BLAS
   */
  const int incx = 1;
  int top = 0;
  int *xi = new int[2 * supNo]();
  //int super_max = 64; //tunig parameter for the max size of supernodes TODO: find it in analysis
  //int col_max = n;
  int *map = new int[n]();
  double *contribs = new double[super_max * col_max]();
  double *trn_diag = new double[super_max * col_max]();
  int info;
  double one[2], zero[2];
  one[0] = 1.0;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
  one[1] = 0.;
  zero[0] = 0.;     /* BETA for *syrk, *herk, and *gemm */
  zero[1] = 0.;

  for (int mn = 0; mn < mod_indices.size(); ++mn) {
   int s = mod_indices[mn] + 1;
   int curCol = s != 0 ? blockSet[s - 1] : 0;
   int nxtCol = blockSet[s];
   int supWdt = nxtCol - curCol;
   int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
   for (int i = Li_ptr[curCol], cnt = 0; i < Li_ptr[nxtCol]; ++i) {
    map[lR[i]] = cnt++;//mapping L rows position to actual row idx
   }

   double *src, *cur = &lValues[lC[curCol]], *cur_d = &D[curCol];//pointing to first element of the current supernode
   //Reseting the current supernode.
   for (int i = 0; i < supWdt; ++i) {
    cur_d[i] = 0;
    for (int j = 0; j < nSupR; ++j) {
     cur[i * nSupR + j] = 0.0;
    }
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
   top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
   assert(top >= 0);
   //std::cout<<"--> "<<s-1<<";;";
   //if(s==2){top =2; xi[top] = 0;}
   for (int i = top; i < supNo; ++i) {
    int lSN = xi[i];
    //std::cout<<lSN<<";";
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


    //  if(ndrow3>0) {

#ifdef OPENBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasConjTrans, ndrow3, ndrow1, supWdts, 1.0, srcL, nSNRCur,
                src, nSNRCur, 0.0, contribs+ndrow1, nSupRs);
#else
    SYM_DGEMM("N","C",&ndrow3,&ndrow1,&supWdts,one,srcL,&nSNRCur,
                        src,&nSNRCur,zero,contribs+ndrow1,&nSupRs );
#endif



//   }

    //copying contrib to L
    for (int i = 0; i < ndrow1; ++i) {//Copy contribs to L
     int col = map[lR[Li_ptr_cSN + i + lb]];//col in the SN
     //double ddiag = 1.0 ;/// D[col];
     for (int j = i; j < nSupRs; ++j) {
      int cRow = lR[Li_ptr_cSN + j + lb];//corresponding row in SN
      //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
      cur[col * nSupR + map[cRow]] -= contribs[i * nSupRs + j];
      //std::cout<<contribs[i*nSupRs+j]<<";";
     }
     // std::cout<<"\n";
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

  return true;
 }
}