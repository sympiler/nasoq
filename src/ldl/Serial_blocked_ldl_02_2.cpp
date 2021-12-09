//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/ldl/Serial_blocked_ldl_02_2.h"

#include <cassert>
#include <chrono>
#include <vector>



#include "nasoq/common/Reach.h"
#include "nasoq/common/Sym_BLAS.h"

namespace nasoq {

 bool ldl_left_sn_02_v2(int n, int *c, int *r, double *values, size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                        double *D, int *blockSet, int supNo, double *timing, int *aTree, int *cT, int *rT, int *col2Sup,
                        int super_max, int col_max, int &nbpivot, int *perm_piv, int *atree_sm, double threshold) {
#if defined(OPENBLAS) && defined(NASOQ_USE_CLAPACK)
  using nasoq::clapacke::LAPACKE_dlapmt;
  using nasoq::clapacke::LAPACKE_dsytrf;
#endif

  /*
   * For timing using BLAS
   */
  std::chrono::time_point<std::chrono::system_clock> start, end1;
  double duration1 = 0;
  const int incx = 1;
  int top = 0;
  int *xi = new int[3 * supNo]();
  int *swap_full = new int[n]();
  std::vector<int> perm_req;
  //int col_max = n;
  int *map = new int[n]();
  double *contribs = new double[super_max * col_max]();
  double *trn_diag = new double[super_max * col_max]();
  int *ws = new int[3 * super_max];
  int *ipiv = new int[super_max]();
  int info;
  double one[2], zero[2];
  one[0] = 1.0;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
  one[1] = 0.;
  zero[0] = 0.;     /* BETA for *syrk, *herk, and *gemm */
  zero[1] = 0.;

  for (int s = 1; s <= supNo; ++s) {

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
     lValues[lC[i] + map[r[j]]] = values[j];
     //std::cout<<i<<";"<<r[j]<<";"<<values[j]<<"\n";
    }
   }
   double *src, *cur = &lValues[lC[curCol]];//pointing to first element of the current supernode
   top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
   assert(top >= 0);
   //int *lbs = xi+supNo, *ubs = xi + 2*supNo;//To use for row permutation
   //if(s==2){top =2; xi[top] = 0;}
   for (int i = top; i < supNo; ++i) {
    int lSN = xi[i];
    int nSupRs = 0;

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
      //lbs[i] = lb;
      sw = false;
     }
     if (lR[j] < curCol + supWdt && !sw) {
      ub = j - Li_ptr_cSN;
      //ubs[i] = ub;
     }
    }
    nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
    int ndrow1 = ub - lb + 1;
    int ndrow3 = nSupRs - ndrow1;
    src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
    double *srcL = &lValues[lC[cSN] + ub + 1];
    blocked_2by2_mult(supWdts, nSupRs, &D[cSN], src, trn_diag, nSNRCur, n);
#ifdef OPENBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasConjTrans, nSupRs, ndrow1, supWdts, 1.0, trn_diag, nSupRs,
                src, nSNRCur, 0.0, contribs, nSupRs);
#else
    SYM_DGEMM("N", "C", &nSupRs, &ndrow1, &supWdts, one, trn_diag, &nSupRs,
          src, &nSNRCur, zero, contribs, &nSupRs);
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
/*     if ( cRow == 78){
      std::cout<<"\n====="<<cSN<<"|| "<< cRow<<";;"<<contribs[i*nSupRs+j]<<";;"
               <<cur[col*nSupR+map[cRow]]<<";;"<<"\n";
     }*/
     }
    }
   }
   LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', supWdt, cur, nSupR, ipiv);
   int is_perm = reorder_after_sytrf(supWdt, cur, nSupR, ipiv,
                                     &perm_piv[curCol], &D[curCol], n, &swap_full[curCol], ws + supWdt);
   // re-order the columns of the super-node
   int rowNo = nSupR - supWdt;
   for (int m = 0; m < supWdt; ++m) {
    perm_piv[curCol + m]++;
   }

   if (is_perm) {
    LAPACKE_dlapmt(LAPACK_COL_MAJOR, 1, rowNo, supWdt, &cur[supWdt], nSupR, &perm_piv[curCol]);
    perm_req.push_back(s);
   }

   //reordering row
   for (int k1 = 0; k1 < supWdt; ++k1) {
    perm_piv[curCol + k1] += (curCol - 1);
    // perm_piv++;
   }
   for (int l = 0; l < supWdt; ++l) {
    D[curCol + l] = cur[l + l * nSupR];
    cur[l + l * nSupR] = 1.0;
   }
#ifdef OPENBLAS
   cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasConjTrans, CblasUnit, rowNo, supWdt, 1.0,
               cur, nSupR, &cur[supWdt], nSupR);
#else
   SYM_DTRSM("R", "L", "C", "U", &rowNo, &supWdt, one,
         cur, &nSupR, &cur[supWdt], &nSupR);
#endif

   blocked_2by2_solver(supWdt, &D[curCol], &cur[supWdt], rowNo, nSupR, n);
  }

/* std::chrono::time_point<std::chrono::system_clock> start, end1;
 start = std::chrono::system_clock::now();*/
  for (int k = 0; k < super_max; ++k) {
   ws[k] = 0;
  }
  row_reordering(supNo, lC, Li_ptr, lR, blockSet, atree_sm, cT, rT, col2Sup,
                 lValues, perm_req, swap_full, xi, map, ws, contribs);
/* end1 = std::chrono::system_clock::now();
 std::chrono::duration<double> elapsed_seconds = end1-start;
 double duration1=elapsed_seconds.count();
 std::cout<<"++ " <<duration1<<"\n";*/
  delete[]contribs;
  delete[]trn_diag;
  delete[]xi;
  delete[]map;
  delete[]ws;
  delete[]swap_full;
  delete []ipiv;
  return true;
 }
}