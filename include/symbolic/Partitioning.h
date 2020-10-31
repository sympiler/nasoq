//
// Created by kazem on 9/25/17.
//

#ifndef CHOLOPENMP_PARTITIONING_H
#define CHOLOPENMP_PARTITIONING_H

#include <vector>

namespace nasoq {
 int partitioning(int inSize, int *inTree, int *inCost,
                  int *inChildPtr, int *inChildNo,//Children list
                  int *nChild, int n, int partitionNum,
   /*Outputs*/
                  int &outSize, int *outCost,
                  int *outNode2Par, std::vector<std::vector<int>> &parList);
}
#endif //CHOLOPENMP_PARTITIONING_H
