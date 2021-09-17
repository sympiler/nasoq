//
// Created by kazem on 1/8/18.
//

#ifndef CHOLOPENMP_INSPECTIONLEVEL_06_H
#define CHOLOPENMP_INSPECTIONLEVEL_06_H

#include <cstddef>

namespace nasoq {
 int getCoarseLevelSet_6(size_t n, const int *eTree, const int *blk2col,
                         int &finaLevelNo, int *&finaLevelPtr, int *&parLevelSet,
                         int &partNo, int *&finalPartPtr, int *&finalNodePtr,
                         int innerParts, int minLevelDist, int divRate,
                         double *nodeCost);
}
#endif //CHOLOPENMP_INSPECTIONLEVEL_06_H
