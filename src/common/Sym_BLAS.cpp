//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/Sym_BLAS.h"

#include <cassert>

//#include "common_interface.h"

namespace nasoq {

 void sym_sytrf(double *A, int n, const int stride, int *nbpivot, double critere) {
  int k;
  double one = 1.0;
#ifdef OPENBLAS
  blasint  iun = 1;
#else
  const int iun = 1;
#endif
  double *tmp, *tmp1;

/* for (int i = 0; i < n; ++i) {
  for (int j = 0; j < n; ++j) {
   std::cout<<A[i*stride+j]<<";";
  }
  std::cout<<"\n";
 }
 for (int j = 0; j < n; ++j) {
  for (int i = 0; i < j; ++i) {
   A[j*stride+i] = 0;
  }
 }*/
  //std::cout<<"diags \n";
  for (k = 0; k < n; k++) {
   tmp = A + k * (stride + 1);
#if 0 //def USE_CSC
   if (std::abs(*tmp)<=critere)
   {
    (*tmp) = critere;
    (*nbpivot)++;
   }
#endif
   tmp1 = tmp + 1; // static PASTIX_FLOAT fun   = 1.0;    static int iun = 1;
   //VEC_SCAL (n-k-1,(one/(*tmp)),tmp1,iun); // Computes the product of a vector by a scalar.
   int tmp_dim = n - k - 1;
   double sca_tmp = double(one / (*tmp));
/*  std::cout<<std::setprecision(48)<<sca_tmp<<";"<<*tmp<<" -- \n";
  for (int i = 0; i < n; ++i) {
   std::cout<<std::setprecision(48)<<tmp1[i]<<";";
  }
  std::cout<<"~~\n";*/
//#ifdef OPENBLAS
   //cblas_dscal(tmp_dim, sca_tmp, tmp1, iun);
//#else
#ifdef OPENBLAS
cblas_dscal(tmp_dim,sca_tmp,tmp1,iun);
#else
   SYM_DSCAL(&tmp_dim, &sca_tmp, tmp1, &iun);
#endif

//#endif
/*  for (int i = 0; i < n; ++i) {
   std::cout<<std::setprecision(48)<<tmp1[i]<<";";
  }
  std::cout<<"++\n";*/
   int dimx = n - k - 1;
   double diag = -(*tmp);
   double *tmp1_stride = tmp1 + stride;
/*  std::cout<<dimx<<":"<<diag<<":"<<*tmp1<<":"<<iun<<":"<<
  *tmp1_stride<<":"<<stride<<" : \n";*/
#ifdef OPENBLAS
   blasint  st = stride;
   cblas_dsyr(CblasColMajor, CblasLower, dimx, diag, tmp1, iun, tmp1_stride, st); //  ?syr Performs a rank-1 update of a symmetric matrix.
  // dsyr_("L", &dimx, &diag, tmp1, &iun, tmp1_stride, &st); //  ?syr Performs a rank-1 update of a symmetric matrix.
#else
   dsyr("L", &dimx, &diag, tmp1, &iun, tmp1_stride, &stride); //  ?syr Performs a rank-1 update of a symmetric matrix.
#endif
/*  for (int i = 0; i < dimx; ++i) {
   std::cout<<tmp1[i]<<";";
  }
  std::cout<<"=== \n";*/
/*  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
    std::cout<<A[i*stride+j]<<";";
   }
   std::cout<<"\n";
  }
  for (int j = 0; j < n; ++j) {
   for (int i = 0; i < j; ++i) {
    A[j*stride+i] = 0;
   }
  }*/
  }
 }

 double dot(int n, double *a, double *b) {
  double result = 0.0;
  for (int i = 0; i < n; ++i) {
   result += (a[i] * b[i]);
  }
  return result;
 }

 void swap_vector(int n, double *a, double *b, int lda) {
  double tmp = 0;
  if (lda == 1) {
   for (int i = 0; i < n; ++i) {// TODO: might need to use a tmp vector for efficincy
    tmp = *(a + i);
    *(a + i) = *(b + i);
    *(b + i) = tmp;
   }
  } else {
   for (int i = 0; i < n; ++i) {
    tmp = *(a + i * lda);
    *(a + i * lda) = *(b + i * lda);
    *(b + i * lda) = tmp;
   }
  }
 }

 void swap_int(int &a, int &b) {
  int tmp = a;
  a = b;
  b = tmp;
 }

 int reorder_after_sytrf(int n, double *a, int lda, int *ipiv, int *perm, double *D, int lda_D, int *swap_vec, int *ws) {
  //Find out actual reordering
  int is_permuted = 0;
  //int *swap_vec = ws;
  int *pivots = ws;
  int *tmpsort = ws + n;
  bool skip_nxt = false;
  for (int j = 0; j < n; ++j) {
   swap_vec[j] = j;//no swap
   tmpsort[j] = j;
  }
  for (int i = 0; i < n; ++i) {
   if (skip_nxt) {
    pivots[i] = 0;
    skip_nxt = false;
    continue;
   }
   int cur_val = ipiv[i]; //1-based indices
   if (cur_val > 0) {
    if (cur_val != i + 1) { // if yes, permutation required
     swap_vec[i] = swap_vec[cur_val - 1];
     is_permuted = 1;
    }
    pivots[i] = 1; //1x1 pivoting
   } else { // 2x2 pivoting
    if (-cur_val != i + 2) {//if yes, perm required
     swap_vec[i + 1] = swap_vec[-cur_val - 1];
     is_permuted = 1;
    }
    pivots[i] = 2; //2x2 pivoting
    D[i + lda_D] = *(a + i * lda + i + 1);
    a[i * lda + i + 1] = 0.0;
    skip_nxt = true;
   }
  }
  // look for 2 indices
/*  for (int l = 0; l < n; ++l) {
   if(pivots[l] == 2){

    D[l+lda_D] = *(a + l*lda + l + 1);
    a [l*lda + l +1] = 0.0;
   }
  }*/
  if (is_permuted) {
   //starting from the end and apply permutation.
   int swap_len = 0;
   for (int k = 0; k < n; ++k) {
    int cur_val = swap_vec[k];
    if (k != cur_val) { //if yes, permutation is needed
     double *col_s1 = a + k; // beginning of row k
     double *col_s2 = a + cur_val; // beginning of row cur_val
     swap_len = k;
     if (pivots[k] == 0) {
      //col_s1-=lda;
      //col_s2-=lda;
      swap_len = k - 1;
     }
     swap_vector(swap_len, col_s1, col_s2, lda);
     swap_int(tmpsort[k], tmpsort[cur_val]);

    }
   }
  }
  // find the sorted indices
  for (int m = 0; m < n; ++m) {
   perm[m] = tmpsort[m];
  }
  return is_permuted;
 }

 void shift(int n_col, int n_row, double *a, int lda, const int inc) {
  if (n_row == 0 || n_col == 0)
   return;
  for (int i = 0; i < n_col; ++i) {
   for (int j = 0; j < n_row; ++j) {
    a[i * lda + inc * j] = a[i * lda + inc * (j + 1)];
   }
  }
 }

 void
 swapping_row_in_blockedL_new(int supWdt, int supNo, const size_t *lC, const size_t *Li_ptr, int *lR, const int curCol,
                              const int *blockSet, double *lValues, int top, int *xi, int *lbs, int *ubs, int *cur_perm,
                              int *cur_swap) {
  double *tmp = new double[supNo]();
  int *map_sup = new int[supWdt]();
  const int minus_one = -1;
  const int one = 1;
// std::cout<<"\nsno: "<<curCol;
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
   //std::cout<<"\nlb,ub "<<lbs[i]<<","<<ubs[i];
   nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
   int num_overlap_rows = ubs[i] - lbs[i];
   double *L_last_row = &lValues[lC[cSN] + ubs[i]]; //beginning of last row
   //double *tmp = &lValues[lC[cNSN-1]]; //FIXME double check
   for (int j = Li_ptr_cSN + lbs[i], cnt = 0; j < Li_ptr_cSN + ubs[i] + 1; ++j, ++cnt) {
    map_sup[lR[j] - curCol] = cnt + 1; //to make sure zero means nothing
/*   if(lR[j] == 207 || lR[j] == 180)
    printf("sss");*/
   }
   for (int m = 0; m < supWdt; ++m) {
    int src_row = m;
    int dst_row = cur_swap[m]; //dst_row always greater than src_row
    if (src_row != dst_row) {
     if (map_sup[src_row] > 0 && map_sup[dst_row] > 0) {//both row exist
      double *row1 = &lValues[lC[cSN] + lbs[i] + map_sup[src_row] - 1];
      double *row2 = &lValues[lC[cSN] + lbs[i] + map_sup[dst_row] - 1];
      swap_vector(supWdts, row1, row2, nSupRs);
      continue;
     } else if (map_sup[src_row] > 0 && map_sup[dst_row] == 0) {
      //finding the closest index to shift
      double *L_src_row = &lValues[lC[cSN] + lbs[i] + map_sup[src_row] - 1];
      for (int kk = 0; kk < supWdts; ++kk) {
       tmp[kk] = L_src_row[nSupRs * kk];
      }
      int first_idx = 0;
      int real_dst_row = curCol + dst_row;
      for (int j = Li_ptr_cSN + lbs[i] + map_sup[src_row] - 1; j < Li_ptr_cSN + ubs[i]; ++j) {
       if (lR[j] > real_dst_row ||
           j == Li_ptr_cNSN - 1) // case of last row
        break;
       first_idx++;
       //shifting lR up as well
       lR[j] = lR[j + 1];
       map_sup[lR[j] - curCol]--;
      }
      // shift up from lR,[srch_row:cur_row_in_blk-1] to lR,[srch_row+1:cur_row_in_blk]
      shift(supWdts, first_idx, L_src_row, nSupRs, one);
      //lR[Li_ptr_cSN + lbs[i] + cnt2] = tmp_row;
      lR[Li_ptr_cSN + lbs[i] + map_sup[src_row] + first_idx - 1] = real_dst_row;
      // copy tmp into L[srch_row]
      if (first_idx > 0) {
       int offset2 = map_sup[src_row] + first_idx - 1;
       double *L_dst_row = &lValues[lC[cSN] + lbs[i] + offset2];
       for (int kk = 0; kk < supWdts; ++kk) {
        L_dst_row[nSupRs * kk] = tmp[kk];
       }
      }
      //tmp[supWdts-1] = lastcol_diag;//diagonal is corrupted, fix it!
      map_sup[dst_row] = map_sup[src_row] + first_idx;
      map_sup[src_row] = 0; // src row is swapped, make it invalid
      continue;
     } else if (map_sup[src_row] == 0 && map_sup[dst_row] > 0) {
      //finding the closest index to shift
      double *L_dst_row = &lValues[lC[cSN] + lbs[i] + map_sup[dst_row] - 1];
      for (int kk = 0; kk < supWdts; ++kk) {
       tmp[kk] = L_dst_row[nSupRs * kk];
      }
      int first_idx = 0;
      int real_src_row = curCol + src_row;
      for (int j = Li_ptr_cSN + lbs[i] + map_sup[dst_row] - 2; //FIXME put it after if
           j > Li_ptr_cSN + lbs[i] - 1; --j) {
       if (lR[j] < real_src_row)
        break;
       first_idx++;
       //shifting lR down as well
       lR[j + 1] = lR[j];
       map_sup[lR[j] - curCol]++;
      }
      // shift down from lR,[srch_row:cur_row_in_blk-1] to lR,[srch_row+1:cur_row_in_blk]
      shift(supWdts, first_idx, L_dst_row, nSupRs, minus_one);
      //lR[Li_ptr_cSN + lbs[i] + cnt2] = tmp_row;
      lR[Li_ptr_cSN + lbs[i] + map_sup[dst_row] - first_idx - 1] = real_src_row;
      // copy tmp into L[srch_row]
      if (first_idx > 0) { // if any row swapped.
       int offset2 = map_sup[dst_row] - first_idx - 1;
       double *L_src_row = &lValues[lC[cSN] + lbs[i] + offset2];
       for (int kk = 0; kk < supWdts; ++kk) {
        L_src_row[nSupRs * kk] = tmp[kk];
       }
      }
      //tmp[supWdts-1] = lastcol_diag;//diagonal is corrupted, fix it!
      map_sup[src_row] = map_sup[dst_row] - first_idx;
      map_sup[dst_row] = 0; // dest row is swapped, make it invalid
      continue;
     } else {//no src and dst rows exists, do nothing
      continue;
     }
    }
   }
   for (int j = 0; j < supWdt; ++j) {
    map_sup[j] = 0; //reset it for nxt supernode
   }
  }
  delete[]tmp;
  delete[]map_sup;
 }

 void
 row_reordering(int supNo, const size_t *lC, const size_t *Li_ptr, int *lR, const int *blockSet, int *aTree, int *cT,
                int *rT, int *col2Sup, double *lValues, std::vector<int> swap_req, int *full_swap, int *xi, int *map,
                int *ws, double *wsf) {
  int *map_sup = ws;
  int *cur_swap;
  const int minus_one = -1;
  const int one = 1;
  double *tmp = wsf;
  for (int t = 0; t < swap_req.size(); ++t) {
   int s = swap_req[t];
   int curCol = s != 0 ? blockSet[s - 1] : 0;
   //std::cout<<"\nsno: "<<curCol;
   int nxtCol = blockSet[s];
   int supWdt = nxtCol - curCol;
   int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
   for (int i = Li_ptr[curCol], cnt = 0; i < Li_ptr[nxtCol]; ++i) {
    map[lR[i]] = cnt++;//mapping L rows position to actual row idx
   }
   double *src, *cur = &lValues[lC[curCol]];//pointing to first element of the current supernode
   int top = ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
   assert(top >= 0);
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
      sw = false;
     }
     if (lR[j] < curCol + supWdt && !sw) {
      ub = j - Li_ptr_cSN;
     }
    }
    //std::cout<<"\nlb,ub "<<lb<<","<<ub;
    nSupRs = Li_ptr_cNSN - Li_ptr_cSN;
    int ndrow1 = ub - lb + 1;
    int ndrow3 = nSupRs - ndrow1;
    src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
    double *srcL = &lValues[lC[cSN] + ub + 1];
    int num_overlap_rows = ub - lb;
    double *L_last_row = &lValues[lC[cSN] + ub]; //beginning of last row
    //double *tmp = &lValues[lC[cNSN-1]]; //FIXME double check
    for (int j = Li_ptr_cSN + lb, cnt = 0; j < Li_ptr_cSN + ub + 1; ++j, ++cnt) {
     map_sup[lR[j] - curCol] = cnt + 1; //to make sure zero means nothing
    }
    cur_swap = &full_swap[curCol];
    for (int m = 0; m < supWdt; ++m) {
     int src_row = m;
     int dst_row = cur_swap[m]; //dst_row always greater than src_row
     if (src_row != dst_row) {
      if (map_sup[src_row] > 0 && map_sup[dst_row] > 0) {//both row exist
       double *row1 = &lValues[lC[cSN] + lb + map_sup[src_row] - 1];
       double *row2 = &lValues[lC[cSN] + lb + map_sup[dst_row] - 1];
       swap_vector(supWdts, row1, row2, nSupRs);
       continue;
      } else if (map_sup[src_row] > 0 && map_sup[dst_row] == 0) {
       //finding the closest index to shift
       double *L_src_row = &lValues[lC[cSN] + lb + map_sup[src_row] - 1];
       for (int kk = 0; kk < supWdts; ++kk) {
        tmp[kk] = L_src_row[nSupRs * kk];
       }
       int first_idx = 0;
       int real_dst_row = curCol + dst_row;
       // see how many nonzero exist after src_row
       for (int j = Li_ptr_cSN + lb + map_sup[src_row] - 1; j < Li_ptr_cSN + ub + 1; ++j) {
        if (lR[j] > real_dst_row)
         break;
        if (j == Li_ptr_cNSN - 1) {// case of last row
         first_idx++;//increase to ensure swap happens
         break;
        }
        first_idx++;
       }
       if (first_idx <= 1) {//no nonzero and no row to swap, just renaming
        assert(lR[Li_ptr_cSN + lb + map_sup[src_row] - 1] == curCol + src_row);
        first_idx = 0;
        lR[Li_ptr_cSN + lb + map_sup[src_row] - 1] = real_dst_row;
        map_sup[dst_row] = map_sup[src_row];
        map_sup[src_row] = 0;
        continue;
       } else {
        first_idx = 0;
        for (int j = Li_ptr_cSN + lb + map_sup[src_row] - 1; j < Li_ptr_cSN + ub + 1; ++j) {
         if (lR[j + 1] > real_dst_row ||
             j == Li_ptr_cNSN - 1) // case of last row
          break;
         first_idx++;
         //shifting lR up as well
         lR[j] = lR[j + 1];
         map_sup[lR[j] - curCol]--;
        }
        // shift up from lR,[srch_row:cur_row_in_blk-1] to lR,[srch_row+1:cur_row_in_blk]
        shift(supWdts, first_idx, L_src_row, nSupRs, one);
        //lR[Li_ptr_cSN + lbs[i] + cnt2] = tmp_row;
        lR[Li_ptr_cSN + lb + map_sup[src_row] + first_idx - 1] = real_dst_row;
        // copy tmp into L[srch_row]
        int offset2 = map_sup[src_row] + first_idx - 1;
        double *L_dst_row = &lValues[lC[cSN] + lb + offset2];
        for (int kk = 0; kk < supWdts; ++kk) {
         L_dst_row[nSupRs * kk] = tmp[kk];
        }

        //tmp[supWdts-1] = lastcol_diag;//diagonal is corrupted, fix it!
        map_sup[dst_row] = map_sup[src_row] + first_idx;
        map_sup[src_row] = 0; // src row is swapped, make it invalid
        continue;
       }
      } else if (map_sup[src_row] == 0 && map_sup[dst_row] > 0) {
       //finding the closest index to shift
       double *L_dst_row = &lValues[lC[cSN] + lb + map_sup[dst_row] - 1];
       for (int kk = 0; kk < supWdts; ++kk) {
        tmp[kk] = L_dst_row[nSupRs * kk];
       }
       int first_idx = 0;
       int real_src_row = curCol + src_row;
       for (int j = Li_ptr_cSN + lb + map_sup[dst_row] - 2; //FIXME put it after if
            j > Li_ptr_cSN + lb - 1; --j) {
        if (lR[j] < real_src_row)
         break;
        first_idx++;
        //shifting lR down as well
        lR[j + 1] = lR[j];
        map_sup[lR[j] - curCol]++;
       }
       // shift down from lR,[srch_row:cur_row_in_blk-1] to lR,[srch_row+1:cur_row_in_blk]
       shift(supWdts, first_idx, L_dst_row, nSupRs, minus_one);
       //lR[Li_ptr_cSN + lbs[i] + cnt2] = tmp_row;
       lR[Li_ptr_cSN + lb + map_sup[dst_row] - first_idx - 1] = real_src_row;
       // copy tmp into L[srch_row]
       if (first_idx > 0) { // if any row swapped.
        int offset2 = map_sup[dst_row] - first_idx - 1;
        double *L_src_row = &lValues[lC[cSN] + lb + offset2];
        for (int kk = 0; kk < supWdts; ++kk) {
         L_src_row[nSupRs * kk] = tmp[kk];
        }
       }
       //tmp[supWdts-1] = lastcol_diag;//diagonal is corrupted, fix it!
       map_sup[src_row] = map_sup[dst_row] - first_idx;
       map_sup[dst_row] = 0; // dest row is swapped, make it invalid
       continue;
      } else {//no src and dst rows exists, do nothing
       continue;
      }
     }
    }
    for (int j = 0; j < supWdt; ++j) {
     map_sup[j] = 0; //reset it for nxt supernode
    }
   }
  }
 }

 void blocked_2by2_solver(int n, double *D, double *rhs, int n_rhs, int lda, int lda_d) {
#ifdef OPENBLAS
  blasint  iun = 1;
#else
  const int iun = 1;
#endif
  for (int i = 0; i < n; ++i) {
   if (D[i + lda_d] == 0) { // simple scaling
    assert(D[i] != 0);
    double tmp = 1.0 / D[i];
#ifdef OPENBLAS
    cblas_dscal(n_rhs, tmp, rhs + i * lda, iun);
#else
    SYM_DSCAL(&n_rhs, &tmp, rhs + i * lda, &iun);
#endif
   } else {//it is 2x2, Cremer rule
    // D[i+n] == D[i+n+1], symm matrix
    //assert(D[i+lda_d]==D[i+lda_d+1]);
    double subdiag = D[i + lda_d];
    double determinant = D[i] * D[i + 1] - subdiag * subdiag;
    double one_over_det = 1.0 / determinant;
    for (int j = 0; j < n_rhs; ++j) {
     double x1 = rhs[i * lda + j];
     double x2 = rhs[(i + 1) * lda + j];
     rhs[i * lda + j] = x1 * D[i + 1] - x2 * subdiag;
     rhs[(i + 1) * lda + j] = x2 * D[i] - x1 * subdiag;
    }
#ifdef OPENBLAS
    cblas_dscal(n_rhs, one_over_det, rhs + i * lda, iun);
    cblas_dscal(n_rhs, one_over_det, rhs + (i + 1) * lda, iun);
#else
    SYM_DSCAL(&n_rhs, &one_over_det, rhs + i * lda, &iun);
    SYM_DSCAL(&n_rhs, &one_over_det, rhs + (i + 1) * lda, &iun);
#endif
    i++;//skip next col since it is part of 2x2 pivoting.
   }
  }
 }

 void blocked_2by2_solver_update(int n, double *D, double *rhs, int n_rhs, int lda, int lda_d, int *mask) {
#ifdef OPENBLAS
  blasint  iun = 1;
#else
  const int iun = 1;
#endif

  for (int i = 0; i < n; ++i) {
   if (mask[i])
    continue;
   if (D[i + lda_d] == 0) { // simple scaling
    assert(D[i] != 0);
    double tmp = 1.0 / D[i];
#ifdef OPENBLAS
    cblas_dscal(n_rhs, tmp, rhs + i * lda, iun);
#else
    SYM_DSCAL(&n_rhs, &tmp, rhs + i * lda, &iun);
#endif
   } else {//it is 2x2, Cremer rule
    // D[i+n] == D[i+n+1], symm matrix
    //assert(D[i+lda_d]==D[i+lda_d+1]);
    double subdiag = D[i + lda_d];
    double determinant = D[i] * D[i + 1] - subdiag * subdiag;
    double one_over_det = 1.0 / determinant;
    for (int j = 0; j < n_rhs; ++j) {
     double x1 = rhs[i * lda + j];
     double x2 = rhs[(i + 1) * lda + j];
     rhs[i * lda + j] = x1 * D[i + 1] - x2 * subdiag;
     rhs[(i + 1) * lda + j] = x2 * D[i] - x1 * subdiag;
    }
#ifdef OPENBLAS
    cblas_dscal(n_rhs, one_over_det, rhs + i * lda, iun);
    cblas_dscal(n_rhs, one_over_det, rhs + (i + 1) * lda, iun);
#else
    SYM_DSCAL(&n_rhs, &one_over_det, rhs + i * lda, &iun);
    SYM_DSCAL(&n_rhs, &one_over_det, rhs + (i + 1) * lda, &iun);
#endif
    i++;//skip next col since it is part of 2x2 pivoting.
   }
  }
 }

 void blocked_2by2_mult(int n, int m, double *D, double *src, double *dst, int lda, int lda_d) {
  int iun = 1;
  for (int l = 0; l < n;) {
   if (D[l + lda_d] == 0) { // simple scaling
    double tmp = D[l];
#if 0
    cblas_dcopy(m,&src[l*lda],iun, &dst[l*m],iun);
    dscal(&m,&tmp,dst+l*m,&iun);
#endif
#if 1 //replaced with blas
    for (int l1 = 0; l1 < m; ++l1) {
     dst[l * m + l1] = tmp * src[l * lda + l1];
    }
#endif
    l++;
   } else {
    double d1 = D[l];
    double d2 = D[l + 1];
    double tmp_d = D[l + lda_d];
    for (int l1 = 0; l1 < m; ++l1) {
     dst[l * m + l1] = d1 * src[l * lda + l1] + tmp_d * src[(l + 1) * lda + l1];
     dst[(l + 1) * m + l1] = tmp_d * src[l * lda + l1] + d2 * src[(l + 1) * lda + l1];
    }
    l += 2;
   }
  }
 }
}