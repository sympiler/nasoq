//
// Created by Kazem on 1/21/19.
//

#ifndef PROJECT_PARALLEL_UPDATE_LDL_02_2_H
#define PROJECT_PARALLEL_UPDATE_LDL_02_2_H

#include <cstddef>

namespace nasoq {
//
// Created by kazem on 12/12/18.
//



#undef TIMING
#undef TLAST
#undef TIMING1
#undef BLASTIMING


 bool update_ldl_left_sn_parallel_02(int n, const int *c, const int *r, const double *values,
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
                                     const int super_max, const int col_max, int &nbpivot,
                                     int *perm_piv,
                                     bool *marked, double threshold = 1e-13);

}
#endif //PROJECT_PARALLEL_UPDATE_LDL_02_2_H
