//
// Created by kazem on 10/8/17.
//

#ifndef CHOLOPENMP_TREEUTILS_H
#define CHOLOPENMP_TREEUTILS_H

#include <vector>
#include <cstdlib>
namespace nasoq {
/*
 * Populates the children
 */
/*int populateChildrenList(int n, int *nChild, int *eTree,
                         int& *childPtr, int& *childNo ){
    auto *childCnt = new int[n]();
    childNo = new int[n];
    childPtr[0]=0;
    for (int k = 1; k < n+1; ++k) {
        childPtr[k] = childPtr[k-1]+nChild[k-1];
    }

    //Populating the children of each node
    for (int l = 0; l < n; ++l) {
        int p = eTree[l];
        if(p>=0){
            childNo[childPtr[p]+childCnt[p]] = l;
            childCnt[p]++;
        }
    }
    delete []childCnt;
    return 0;
}*/

 void populateChildren(int n, const int *eTree, int *childPtr,
                       int *childNo, int *nChild);

/*
 *
 */
 int getNodeDepth(int node, int n, const int *tree, int *weight = NULL);

/*
 *
 */
 int getTreeHeightBruteForce(int n, const int *tree, int *weight = NULL);

 int getTreeHeight(int n, const int *tree, int *nChild1, int *weight = NULL);

/*
 *
 */
 double *computeSubtreeCost(int n, const int *tree, int *nChild, double *weight);


 int getLevelSet(size_t n, const int *inTree, int *levelPtr, int *levelSet);


 int mergeInnerPart(std::vector<std::vector<int>> newLeveledParList,
                    int *inCost, std::vector<std::vector<int>> &mergedLeveledParList,
                    int *outCost,
                    int costThreshold);

 struct subTree {
  double cost;
  std::vector<int> nodeList;
 };

 bool cmpCost(subTree a, subTree b);

 int findMin(double *cost, int size);

 int worstFitBinPack(const std::vector<std::vector<int>> &newLeveledParList,
                     double *inCost,
                     std::vector<std::vector<int>> &mergedLeveledParList,
                     double *outCost, int costThreshold, int numOfBins);

 int heightPartitioning_DAG_Trng(int levelNo,
                                 int *levelPtr,
                                 int *node2Level,
                                 int originalHeight,
                                 int innerParts,
                                 int minLevelDist,
                                 int divRate,
                                 std::vector<int> &innerPartsSize,
                                 std::vector<std::vector<int>> &slackGroups,
                                 double *subTreeCost,
                                 int *partition2Level,
                                 bool sw = false
 );


 int heightPartitioning(int levelNo,
                        int *levelPtr,
                        int *node2Level,
                        int originalHeight,
                        int innerParts,
                        int minLevelDist,
                        int divRate,
                        std::vector<int> &innerPartsSize,
                        std::vector<std::vector<int>> &slackGroups,
                        double *subTreeCost,
                        int *partition2Level,
                        bool sw = false
 );

 void makeSlackedLevelSet(int n, int clusterCnt,
                          int *partition2Level,
                          int originalHeight,
                          std::vector<std::vector<int>> &slackGroups,//out
                          std::vector<std::vector<int>> &slackedLevelSet,//out
                          int *node2Level);
}
#endif //CHOLOPENMP_TREEUTILS_H
