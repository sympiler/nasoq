//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/matrixMatrix/spmm.h"

namespace nasoq {

 void spmdm(int n, size_t *Ap, int *Ai, double *valA, double *D, double *valAD) {
  for (int i = 0; i < n; ++i) {
   if (D[i + n] == 0) {
    double diag = D[i];
    double *ccA = &valA[Ap[i]];
    double *ccAD = &valAD[Ap[i]];
    for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
     *(ccAD++) = diag * *(ccA++);
    }
   } else {
    double d1 = D[i];
    double d2 = D[i + 1];
    double tmp_d = D[i + n];
    valAD[Ap[i]] = d1; //FIXME: A should have a nonzero zero in diagonal for 2x2
    for (int j = Ap[i] + 1, k = Ap[i + 1]; j < Ap[i + 1]; ++j, ++k) {
     valAD[j] = d1 * valA[j] + tmp_d * valA[k];
     valAD[k] = tmp_d * valA[j] + d1 * valA[k];
    }
    i++;
   }
  }
 }
}