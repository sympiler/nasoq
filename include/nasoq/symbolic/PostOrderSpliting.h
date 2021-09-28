//
// Created by kazem on 1/25/18.
//

#ifndef CHOLOPENMP_POSTORDERSPLITING_H
#define CHOLOPENMP_POSTORDERSPLITING_H

#include <vector>

namespace nasoq {
 int postOrderSpliting(int inSize, int *inTree, double *inCost,
                       int *inChildPtr, int *inChildNo,//Children list
                       int *nChild, int n, int partitionNum,
   /*Outputs*/
                       int &outSize, double *outCost,
                       int *outNode2Par, std::vector<std::vector<int>> &parList);
}
#endif //CHOLOPENMP_POSTORDERSPLITING_H
