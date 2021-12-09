//
// Created by kazem on 11/20/18.
//
#ifndef PROJECT_LDLT_CHECK_H
#define PROJECT_LDLT_CHECK_H

#include <cstddef>

#include "nasoq/matrixVector/spmv_CSC.h"

namespace nasoq {
/*
 * Assuming A' is passed which is basically CSR of A
 */
 bool ldlt_check(int n, size_t *Lp, int *Li, double *valL,
                 double *D, size_t *Ap, int *Ai, double *valA);
}
#endif //PROJECT_LDLT_CHECK_H
