//
// Created by Kazem on 12/9/18.
//

#ifndef PROJECT_SERIAL_BLOCKED_LDL01_H
#define PROJECT_SERIAL_BLOCKED_LDL01_H
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include "mkl.h"
#include "Reach.h"
#include "Sym_BLAS.h"

namespace nasoq{

/*
 * LDLT with bunch-kaufman pivoting in each supernode
 * Does row reordering during factorization
 */

bool ldl_left_sn_02(int n, int* c, int* r, double* values,
                     size_t *lC, int * lR, size_t * Li_ptr, double* lValues,
                     double *D,
                     int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                     int *aTree, int *cT, int *rT, int *col2Sup,
#else
  int *prunePtr, int *pruneSet,
#endif
                     int super_max, int col_max, int &nbpivot, int *perm_piv,
                     double threshold=1e-13) {
  /*
   * For timing using BLAS
   */
  const int incx = 1;
  int top = 0;
  int *xi = new int[3 * supNo]();
  //int super_max = 64; //tunig parameter for the max size of supernodes TODO: find it in analysis
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

/*  if(curCol == 36)
   printf("dfgdf");*/
#ifndef PRUNE
   top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
   assert(top >= 0);
   int *lbs = xi + supNo, *ubs = xi + 2 * supNo;//To use for row permutation
   //if(s==2){top =2; xi[top] = 0;}
   for (int i = top; i < supNo; ++i) {
    int lSN = xi[i];
#else
    for (int i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
       int lSN = pruneSet[i];
#endif
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
      lbs[i] = lb;
      sw = false;
     }
     if (lR[j] < curCol + supWdt && !sw) {
      ub = j - Li_ptr_cSN;
      ubs[i] = ub;
     }
    }
    nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
    int ndrow1 = ub - lb + 1;
    int ndrow3 = nSupRs - ndrow1;
    src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
    double *srcL = &lValues[lC[cSN] + ub + 1];
    // multiplying L * D
/*   for (int l = 0; l < supWdts; ++l) {
    double tmp = D[cSN+l];
    for (int l1 = 0; l1 < nSupRs; ++l1) {
     trn_diag[l*nSupRs+l1] = tmp * src[l*nSNRCur+l1];
    }
   }*/
/*if(cSN==8)
 printf("ggggg");*/
    blocked_2by2_mult(supWdts, nSupRs, &D[cSN], src, trn_diag, nSNRCur, n);
/*   if (s == supNo){
    std::cout<<"Col: "<<cSN<<"\n";
    for (int l = 0; l < supWdts; ++l) {
     for (int l1 = 0; l1 < nSupRs; ++l1) {
      std::cout<<trn_diag[l*nSupRs+l1]<<";" ;
     }
     std::cout<<"\n";
    }
   }*/


    //  if(ndrow3>0) {
    /*dgemm("N", "C", &ndrow3, &ndrow1, &supWdts, one, srcL, &nSNRCur,
          src, &nSNRCur, zero, &contribs[ndrow1], &nSupRs);*/
    dgemm("N", "C", &nSupRs, &ndrow1, &supWdts, one, trn_diag, &nSupRs,
          src, &nSNRCur, zero, contribs, &nSupRs);

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

   //dpotrf("L",&supWdt,cur,&nSupR,&info);

/*   for (int m = 0; m < supWdt; ++m) {
    for (int i = m; i < supWdt; ++i) {
     contribs[m*su]
    }
   }

   dspff*/
/*  if (s == supNo) {
   std::cout << "+\n";
   for (int m = 0; m < supWdt; ++m) {
    for (int i = 0; i < supWdt; ++i) {
     std::cout << cur[i + m * nSupR] << ";";
    }
    std::cout << "+\n";
   }
  }*/
   //sym_sytrf ( cur, supWdt, nSupR, &nbpivot, threshold);
   LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', supWdt, cur, nSupR, ipiv);


/*  std::cout<<"+\n";
  for (int m = 0; m < supWdt; ++m) {
   for (int i = 0; i < supWdt; ++i) {
    std::cout<<cur[i+m*nSupR]<<";";
   }
   std::cout<<"+\n";
  }*/
   int is_perm = reorder_after_sytrf(supWdt, cur, nSupR, ipiv, &perm_piv[curCol], &D[curCol], n, ws, ws + supWdt);

/* if (curCol == 77){
  std::cout<<"+\n";
 for (int m = 0; m < supWdt; ++m) {
   for (int i = 0; i < supWdt; ++i) {
    std::cout<<cur[i+m*nSupR]<<";";
   }
   std::cout<<"\n";
  }
  std::cout<<"+\n";
  }*/
#if 0
   std::cout<<"** diag "<< curCol <<"\n";
   for (int m = 0; m < supWdt; ++m) {
    for (int i = 0; i < supWdt; ++i) {
     std::cout<<cur[i+m*nSupR]<<";";
    }
    std::cout<<"\n";
   }
#endif
   //std::cout<<"++ ord "<< curCol <<"\n";
/*  for (int i1 = 0; i1 < supWdt; ++i1) {
   if (ipiv[i1] < 0)
    std::cout<<"\n"<<ipiv[i1]<<";====";
  }*/
   //std::cout<<"\n";

   // re-order the columns of the super-node

   int rowNo = nSupR - supWdt;
   for (int m = 0; m < supWdt; ++m) {
    perm_piv[curCol + m]++;
   }
/*if (curCol == 77){
 std::cout<<"***\n";
 for (int i1 = 0; i1 < supWdt; ++i1) {
  std::cout<<ws[i1]<<";";
 }
 std::cout<<"***\n";
 std::cout<<"***\n";
 for (int k1 = 0; k1 < supWdt; ++k1) {
  std::cout<< perm_piv[curCol+k1]<<";";
 }
 std::cout<<"***\n";
#if 0
 std::cout<<"== offdiag "<< curCol <<"\n";
 if(rowNo>0)
  for (int k = 0; k < supWdt; ++k) {
   for (int i = 0; i < rowNo; ++i) {
    std::cout<<cur[supWdt + k*nSupR + i]<<",";
   }
   std::cout<<"\n";
  }
#endif
}*/

   LAPACKE_dlapmt(LAPACK_COL_MAJOR, 1, rowNo, supWdt, &cur[supWdt], nSupR, &perm_piv[curCol]);
   //reordering row
   for (int k1 = 0; k1 < supWdt; ++k1) {
    perm_piv[curCol + k1] += (curCol - 1);
    // perm_piv++;
   }
/*  if(curCol==77){
   std::cout<<"== offdiag "<< curCol <<"\n";
   if(rowNo>0)
    for (int k = 0; k < supWdt; ++k) {
     for (int i = 0; i < rowNo; ++i) {
      std::cout<<cur[supWdt + k*nSupR + i]<<",";
     }
     std::cout<<"\n";
    }
  }*/
   if (is_perm) {
    /*swapping_row_in_blockedL(supNo, lC, Li_ptr, lR,
     curCol, blockSet, lValues, top, xi, lbs, ubs, perm_piv, ws);*/
    swapping_row_in_blockedL_new(supWdt, supNo, lC, Li_ptr, lR,
                                 curCol, blockSet, lValues,
                                 top, xi, lbs, ubs, perm_piv, ws);
   }
   // Making L*D
/*  for (int l = 0; l < supWdt; ++l) {
   double tmp = cur[l + l*nSupR];
   D[curCol+l] = tmp;
   double *stCol = trn_diag + l*supWdt + l;
   double *curCol = cur + l*nSupR + l;
   *stCol = tmp;
   for (int l1 = 0; l1 < supWdt-l-1; ++l1) {
    *(++stCol) = tmp * *(++curCol);
   }
  }
  dtrsm("R", "L", "C", "N", &rowNo, &supWdt,one,
        trn_diag,&supWdt,&cur[supWdt],&nSupR);

  for (int k = 0; k < supWdt; ++k) {
   cur[k*nSupR+k] = 1.0;
  }*/
   for (int l = 0; l < supWdt; ++l) {
    D[curCol + l] = cur[l + l * nSupR];
    cur[l + l * nSupR] = 1.0;
   }
#if 0
   if(curCol==36)
      if(true){
         std::cout<<"\n+++++\n";
       for (int m = 0; m < supWdt; ++m) {
        for (int i = 0; i < supWdt; ++i) {
         std::cout<<cur[i+m*nSupR]<<";";
        }
        std::cout<<"+\n";
       }
      std::cout<<"== offdiag "<< curCol <<"\n";
     if(rowNo>0)
     for (int k = 0; k < supWdt; ++k) {
      for (int i = 0; i < rowNo; ++i) {
       std::cout<<cur[supWdt + k*nSupR + i]<<",";
      }
      std::cout<<"\n";
     }}
#endif
   dtrsm("R", "L", "C", "U", &rowNo, &supWdt, one,
         cur, &nSupR, &cur[supWdt], &nSupR);
   blocked_2by2_solver(supWdt, &D[curCol], &cur[supWdt], rowNo, nSupR, n);

   //copying 1/Di into D
   /*for (int l = 0; l < supWdt; ++l) {
    D[curCol+l] = one[0] / cur[l + l*nSupR];
   }*/

#if 0
   if(curCol == 36){
   std::cout<<"== offdiag "<< curCol <<"\n";
   if(rowNo>0)
   for (int k = 0; k < supWdt; ++k) {
    for (int i = 0; i < rowNo; ++i) {
     std::cout<<cur[supWdt + k*nSupR + i]<<",";
    }
    std::cout<<"\n";
   }}
#endif
  }

  delete[]contribs;
  delete[]trn_diag;
  delete[]xi;
  delete[]map;
  delete[]ws;

  return true;
 }


}

#endif //PROJECT_SERIAL_BLOCKED_LDL01_H
