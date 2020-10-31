//
// Created by kazem on 11/20/18.
//

#ifndef PROJECT_SPMM_H
#define PROJECT_SPMM_H

#include <cstddef>

namespace nasoq {
/*
 * Sparse matrix times diagonal matrix
 * the pointers can be reused for the output matrix.
 * only new values are returned.
 */
 void spmdm(int n, size_t *Ap, int *Ai, double *valA,
            double *D, double *valAD);

}
#endif //PROJECT_SPMM_H
