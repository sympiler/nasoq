//
// Created by kazem on 11/5/18.
//

#ifndef PROJECT_SERIAL_BLOCKED_LDL_H
#define PROJECT_SERIAL_BLOCKED_LDL_H

#include <cstddef>

namespace nasoq{

/*
 * LDLT factorization without pivoting or static pivoting
 */


bool ldl_left_sn_01 (int n, int* c, int* r, double* values,
                         size_t *lC, int * lR, size_t * Li_ptr, double* lValues,
                         double *D,
                         int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                    int *aTree, int *cT, int *rT, int *col2Sup,
#else
  int *prunePtr, int *pruneSet,
#endif
                    int super_max, int col_max, int &nbpivot,
                     double threshold=1e-13);

}
#endif //PROJECT_SERIAL_BLOCKED_LDL_H
