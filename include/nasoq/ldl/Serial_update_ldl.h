//
// Created by kazem on 1/10/19.
//

#ifndef PROJECT_SERIAL_UPDATE_LDL_H
#define PROJECT_SERIAL_UPDATE_LDL_H

#include <cstddef>
#include <vector>

namespace nasoq {
/*
 * Update-downdate algorithm for LDLT with
 * bunch-kaufman pivoting in each supernode
 * Does row reordering after factorization
 *
 * Int workspace: 4*supermax + 2*n + 3*supNo
 * double workspace: 2 * super_max*col_max
 */

 bool update_ldl_left_sn_02_v2(int n, int *c, int *r, double *values,
                               size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                               double *D,
                               int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                               int *aTree, int *cT, int *rT, int *col2Sup,
#else
   int *prunePtr, int *pruneSet,
#endif
                               std::vector<int> mod_indices,
                               int super_max, int col_max, int &nbpivot, int *perm_piv,
                               int *atree_sm, int *ws_int = NULL, double *ws_dbl = NULL,
                               double threshold = 1e-13);
}
#endif //PROJECT_SERIAL_UPDATE_LDL_H
