//
// Created by Kazem on 4/9/19.
//

#ifndef PROJECT_SERIAL_SIMPLICIAL_LDL_H
#define PROJECT_SERIAL_SIMPLICIAL_LDL_H


namespace nasoq{

int ldl_left_simplicial_01(int n, int* c, int* r, double* values,
                           int* cT, int* rT,
                           int* lC, int* lR, double* &lValues,
                           double *d,
#if 0
  int *prunePtr, int *pruneSet,
#endif
                           int *eTree,
                           double *ws, int *ws_int);


 int ldl_left_simplicial_02(int n, int *c, int *r, double *values,
                            int *cT, int *rT,
                            int *lC, int *lR, double *&lValues,
                            double *d,
#if 0
   int *prunePtr, int *pruneSet,
#endif
                            int *eTree,
                            double *ws, int *ws_int);
}
#endif //PROJECT_SERIAL_SIMPLICIAL_LDL_H
