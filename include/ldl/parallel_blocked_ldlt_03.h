//
// Created by Kazem on 12/31/18.
//

#ifndef PROJECT_SERIAL_BLOCKED_LDL_03_H
#define PROJECT_SERIAL_BLOCKED_LDL_03_H

#include <cstddef>

namespace nasoq {
#undef TIMING
#undef TLAST


 bool ldl_left_sn_parallel_03(int n, const int *c, const int *r, const double *values,
                              const size_t *lC, int *lR, const size_t *Li_ptr, double *lValues,
                              double *D,
                              const int *blockSet, const int supNo, double *timing,
#ifndef PRUNE
                              int *aTree, int *cT, int *rT, int *col2Sup,
#else
   int *prunePtr, int *pruneSet,
#endif
                              const int nLevels, const int *levelPtr, const int *levelSet,
                              const int nPar, const int *parPtr, const int *partition,
                              const int chunk, const int threads,
                              const int super_max, const int col_max, int &nbpivot, int *perm_piv,
                              double threshold = 1e-13);
}
#endif //PROJECT_SERIAL_BLOCKED_LDL_03_H
