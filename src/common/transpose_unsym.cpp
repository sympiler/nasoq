//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/transpose_unsym.h"

#include <algorithm>

namespace nasoq {

 int
 transpose_unsym(size_t row, size_t col, int *Ap, int *Ai, double *Ax, size_t &rowT, size_t &colT, int *&ATp, int *&ATi,
                 double *&ATx) {
  rowT = col;
  colT = row;
  if (row == 0) {
   ATp = NULL;
   ATi = NULL;
   ATx = NULL;
   return 1;
  }
  int *col_cnt = new int[row]();//cols will be rows
  for (int i = 0; i < Ap[col]; ++i) {
   col_cnt[Ai[i]]++;
  }
  ATp = new int[colT + 1];
  ATp[0] = 0;
  //Determining column pointer of AT
  for (int j = 1; j < colT + 1; ++j) {
   ATp[j] = ATp[j - 1] + col_cnt[j - 1];
  }
  std::fill_n(col_cnt, colT, 0);
  //Determining row pointer of AT
  ATi = new int[ATp[colT]];
  ATx = new double[ATp[colT]];
  for (int k = 0; k < col; ++k) {
   for (int i = Ap[k]; i < Ap[k + 1]; ++i) {
    int beg = Ai[i];
    ATi[ATp[beg] + col_cnt[beg]] = k;
    ATx[ATp[beg] + col_cnt[beg]] = Ax[i];
    col_cnt[beg]++;
   }
  }
  delete[]col_cnt;
  return 1;
 }

}
