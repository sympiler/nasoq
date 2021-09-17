//
// Created by kazem on 3/13/19.
//

#ifndef PROJECT_TRANSPOSE_UNSYM_H
#define PROJECT_TRANSPOSE_UNSYM_H

#include <cstddef>

namespace nasoq {
 int transpose_unsym(size_t row, size_t col, int *Ap, int *Ai, double *Ax,
                     size_t &rowT, size_t &colT, int *&ATp, int *&ATi,
                     double *&ATx);
}
#endif //PROJECT_TRANSPOSE_UNSYM_H
