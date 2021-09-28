//
// Created by kazem on 4/12/19.
//

#ifndef PROJECT_PARALLEL_SIMPLICIAL_LDL_H
#define PROJECT_PARALLEL_SIMPLICIAL_LDL_H

namespace nasoq{
int ldl_parallel_left_simplicial_01 (int n, int* c, int* r, double* values,
                           int* cT, int* rT,
                           int* lC, int* lR, double* &lValues,
                           double *d,
#if 0
  int *prunePtr, int *pruneSet,
#endif
                           int *eTree,
                           int nLevels, int *levelPtr,
                           int nPar, int *parPtr, int *partition);
}
#endif //PROJECT_PARALLEL_SIMPLICIAL_LDL_H
