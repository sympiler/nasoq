//
// Created by kazem on 4/11/19.
//

#ifndef PROJECT_SERIAL_UPDATE_SIMPLICIAL_LDL_H
#define PROJECT_SERIAL_UPDATE_SIMPLICIAL_LDL_H

#include <vector>

namespace nasoq {

 int update_ldl_left_simplicial_01(int n, int *c, int *r, double *values,
                                   int *cT, int *rT,
                                   int *lC, int *lR, double *&lValues,
                                   double *d,
#if 0
   int *prunePtr, int *pruneSet,
#endif
                                   int *eTree, std::vector<int> mod_indices,
                                   double *ws, int *ws_int);
}

#endif //PROJECT_SERIAL_UPDATE_SIMPLICIAL_LDL_H
