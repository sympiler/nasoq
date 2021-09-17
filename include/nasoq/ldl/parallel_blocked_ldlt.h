//
// Created by kazem on 11/30/18.
//

#ifndef PROJECT_PARALLEL_BLOCKED_LDLT_H
#define PROJECT_PARALLEL_BLOCKED_LDLT_H

#include <cstddef>

#undef TIMING
#undef TLAST
#undef TIMING1
namespace nasoq {
 bool ldl_left_sn_parallel_01(int n, int *c, int *r, double *values,
                              size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                              double *D,
                              int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                              int *aTree, int *cT, int *rT, int *col2Sup,
#else
   int *prunePtr, int *pruneSet,
#endif
                              int nLevels, int *levelPtr, int *levelSet,
                              int nPar, int *parPtr, int *partition,
                              int chunk, int threads,
                              int super_max, int col_max, int &nbpivot,
                              double threshold = 1e-13);

#undef TIMING1


#undef BLASTIMING
}
#endif //PROJECT_PARALLEL_BLOCKED_LDLT_H
