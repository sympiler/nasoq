//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/TreeUtils.h"

#include <algorithm>

#include "nasoq/common/def.h"

namespace nasoq {

 void populateChildren(int n, const int *eTree, int *childPtr, int *childNo, int *nChild) {
  int *childCnt = new int[n]();
  for (int k = 0; k < n; ++k) {
   if (eTree[k] >= 0)
    nChild[eTree[k]]++;
  }
  childPtr[0] = 0;
  for (int k = 1; k < n + 1; ++k) {
   childPtr[k] = childPtr[k - 1] + nChild[k - 1];
  }
  for (int l = 0; l < n; ++l) {
   int p = eTree[l];
   if (p >= 0) {
    childNo[childPtr[p] + childCnt[p]] = l;
    childCnt[p]++;
   }
  }
  delete[]childCnt;
 }

 int getNodeDepth(int node, int n, const int *tree, int *weight) {
  int level = 0;
  if (weight != NULL)
   level += weight[node];
  while (tree[node] >= 0) {
   node = tree[node];
   if (weight == NULL)
    level++;
   else
    level += weight[node];
  }
  //return level+1;
  return level;
 }

 int getTreeHeightBruteForce(int n, const int *tree, int *weight) {
  int maxLen = 0;
  for (int i = 0; i < n; ++i) {
   int ltmp = getNodeDepth(i, n, tree, weight);
   if (ltmp > maxLen)
    maxLen = ltmp;
  }

  return maxLen;
 }

 int getTreeHeight(int n, const int *tree, int *nChild1, int *weight) {
  int maxLen = 0;

  for (int i = 0; i < n; ++i) {
   if (nChild1[i] == 0) {
    int ltmp = getNodeDepth(i, n, tree, weight);
    if (ltmp > maxLen)
     maxLen = ltmp;
   }
  }
  return maxLen;
 }

 double *computeSubtreeCost(int n, const int *tree, int *nChild, double *weight) {
  double *subTreeCost = new double[n];
  for (int i = 0; i < n; ++i) {//each subtree has a node
   subTreeCost[i] = weight[i];
  }
  for (int i = 0; i < n; ++i) {
   int par = tree[i];
   if (par > 0) {
    subTreeCost[par] += subTreeCost[i];
   }
  }
  return subTreeCost;
 }

 int getLevelSet(size_t n, const int *inTree, int *levelPtr, int *levelSet) {
//Naive code for generating level levelSet from ETree, it is not part of
// the inspector. In a real code, we return a node or supernode at
// a time so, it will in O(n)
  int begin = 0, end = (int) n - 1;
  auto *nChild = new int[n]();
  auto *visited = new bool[n]();
  //Counting the number of children
  for (int k = 0; k < n; ++k) {
   if (inTree[k] >= 0)
    nChild[inTree[k]]++;
  }
  int curLevel = 0, curLevelCnt = 0;
  levelPtr[0] = 0;

  while (begin <= end) {
   for (int i = begin; i <= end; ++i) {//For level curLevel
    if (nChild[i] == 0 && !visited[i]) {//if no incoming edge
     visited[i] = true;
     levelSet[curLevelCnt] = i; //add it to current level
     curLevelCnt++;//Adding to level-set
    }
   }
   curLevel++;//all nodes with zero indegree are processed.
   levelPtr[curLevel] = curLevelCnt;
   if (curLevelCnt == n)
    break;
   while (nChild[begin] == 0)
    begin++;
   while (nChild[end] == 0 && begin <= end)
    end--;
   //Updating degrees after removing the nodes
   for (int l = levelPtr[curLevel - 1]; l < levelPtr[curLevel]; ++l) {
    int cc = levelSet[l];
    if (inTree[cc] >= 0)
     nChild[inTree[cc]]--;
   }
  }
#if DEBUG >= 2
  std::cout<<"FinalSet\n";
  for (int l = 0; l < curLevel; ++l) {
   for (int lp = levelPtr[l]; lp < levelPtr[l+1]; ++lp) {
    std::cout<<levelSet[lp]<<",";
   }
   std::cout<<"\n\n";
  }
  std::cout<<"\n";
#endif
  delete[]nChild;
  delete[]visited;
  return curLevel;

 }

 int mergeInnerPart(std::vector<std::vector<int>> newLeveledParList, int *inCost,
                    std::vector<std::vector<int>> &mergedLeveledParList, int *outCost, int costThreshold) {
  int lClusterCnt = 0;
  int partNo = newLeveledParList.size();
#ifdef DEBUG
  for (int j = 0; j < partNo; ++j) {
         std::cout<<inCost[j]<<",";
     }
     std::cout<<"\n";
#endif
  for (int i = 0; i < partNo;) {
   int curCost = 0;
   while (curCost < (1 * costThreshold) && i < partNo) {
    curCost += inCost[i];
    mergedLeveledParList[lClusterCnt].insert(
      mergedLeveledParList[lClusterCnt].end(),
      newLeveledParList[i].begin(), newLeveledParList[i].end());
    i++;
   }
   lClusterCnt++;
   outCost[lClusterCnt] = curCost;
  }
  return lClusterCnt;
 }

 bool cmpCost(subTree a, subTree b) {
  return a.cost > b.cost;
 }

 int findMin(double *cost, int size) {
  double min = INT_MAX;
  int minBin = 0;
  for (int i = 0; i < size; ++i) {
   if (cost[i] < min) {
    min = cost[i];
    minBin = i;
   }
  }
  return minBin;
 }

 int worstFitBinPack(const std::vector<std::vector<int>> &newLeveledParList, double *inCost,
                     std::vector<std::vector<int>> &mergedLeveledParList, double *outCost, int costThreshold,
                     int numOfBins) {

  int lClusterCnt = 0;
  int partNo = newLeveledParList.size();
  //Sorting the subtree list
  std::vector<subTree> partList(partNo);
  for (int i = 0; i < partNo; ++i) {
   partList[i].cost = inCost[i];
   partList[i].nodeList.insert(partList[i].nodeList.begin(),
                               newLeveledParList[i].begin(), newLeveledParList[i].end());
  }
  std::sort(partList.begin(), partList.end(), cmpCost);

#if 0//def DEBUG
  for (int j = 0; j < partNo; ++j) {
    std::cout<<partList[j].cost<<",";
  }
  std::cout<<"\n";
#endif
  int minBin = 0;
  for (int i = 0; i < partNo; i++) {
   minBin = findMin(outCost, numOfBins);
   outCost[minBin] += partList[i].cost;
   mergedLeveledParList[minBin].insert(
     mergedLeveledParList[minBin].end(),
     partList[i].nodeList.begin(), partList[i].nodeList.end());
  }
#if SHOWCOST//def DEBUG
  for (int j = 0; j < numOfBins; ++j) {
   std::cout<<outCost[j]<<";";
  }
  //std::cout<<"\n";
#endif

  return numOfBins;
 }

 int heightPartitioning_DAG_Trng(int levelNo, int *levelPtr, int *node2Level, int originalHeight, int innerParts,
                                 int minLevelDist, int divRate, std::vector<int> &innerPartsSize,
                                 std::vector<std::vector<int>> &slackGroups, double *subTreeCost, int *partition2Level,
                                 bool sw) {
  int levelCut = 0, preLevelCut = 0, finalSeqNodes = 3;
  int lClusterCnt = 0, innerPartsTmp = innerParts;
  //auto partition2Level = new int[levelNo+1]();
  int *accuSlackGroups = new int[levelNo];
  if (levelNo <= minLevelDist) {
   partition2Level[0] = 0;
   partition2Level[1] = levelNo;
   innerPartsSize.push_back(1);
   lClusterCnt = 1;
   return lClusterCnt;
  }
  //assigne the nodes in normal level set
  for (int i = 0; i < levelNo; ++i) {
   /*accuSlackGroups[i] = levelPtr[i+1]-levelPtr[i] +
     slackGroups[i].size();*/
   accuSlackGroups[i] = levelPtr[i + 1] - levelPtr[i];
   assert(accuSlackGroups[i] >= 0);
  }


  partition2Level[0] = 0;


  /*if(minLevelDist<=0)
   minLevelDist=2;//default parameter*/
  int tmp = 0;
  tmp = minLevelDist;
  if (tmp > partition2Level[lClusterCnt] && tmp < levelNo) {
   //Due to tuning parameter we need this check
   partition2Level[++lClusterCnt] = tmp;
   innerPartsTmp = accuSlackGroups[tmp - 1] / 2 > 1 ? accuSlackGroups[tmp - 1] : 1;
   innerPartsSize.push_back(innerPartsTmp);
  }
  tmp += divRate;
  while (tmp < originalHeight - 1) {
   //Ensures a certain number of level in each partition
   innerPartsTmp = accuSlackGroups[tmp - 1];
   innerPartsSize.push_back(innerPartsTmp > 1 ? innerPartsTmp : 1);
   partition2Level[++lClusterCnt] = tmp;
   tmp += divRate;
  }
  partition2Level[++lClusterCnt] = originalHeight + 1;
  innerPartsSize.push_back(1);//The last partition has one element
  delete[]accuSlackGroups;
#if 0
  for (int i1 = 0; i1 <= lClusterCnt; ++i1) {
    std::cout << partition2Level[i1] << ",";
   }
   for (int i1 = lClusterCnt + 1; i1 < 7; ++i1) {
    std::cout << partition2Level[i1] << ",";
   }
  //std::cout<<"\n";
#endif
  return lClusterCnt;
 }

 int
 heightPartitioning(int levelNo, int *levelPtr, int *node2Level, int originalHeight, int innerParts, int minLevelDist,
                    int divRate, std::vector<int> &innerPartsSize, std::vector<std::vector<int>> &slackGroups,
                    double *subTreeCost, int *partition2Level, bool sw) {
  int levelCut = 0, preLevelCut = 0, finalSeqNodes = 3;
  int lClusterCnt = 0, innerPartsTmp = innerParts;
  //auto partition2Level = new int[levelNo+1]();
  if (levelNo <= 2) {
   partition2Level[0] = 0;
   partition2Level[1] = levelNo;
   innerPartsSize.push_back(1);
   lClusterCnt = 1;
   return lClusterCnt;
  }
  int *accuSlackGroups = new int[levelNo];
  //assigne the nodes in normal level set
  for (int i = 0; i < levelNo; ++i) {
   /*accuSlackGroups[i] = levelPtr[i+1]-levelPtr[i] +
     slackGroups[i].size();*/
   accuSlackGroups[i] = levelPtr[i + 1] - levelPtr[i];
   assert(accuSlackGroups[i] >= 0);
  }


  partition2Level[0] = 0;

  if (sw) {
   for (int i = minLevelDist; i < levelNo; i += minLevelDist) {
    lClusterCnt++;
    partition2Level[lClusterCnt] = i;
    //for leaves
    if (accuSlackGroups[i] >= divRate * innerPartsTmp) {
     innerPartsSize.push_back(innerPartsTmp);
    } else { //otherwise a divisor of divRate
     int tmp = accuSlackGroups[i] / 2;
     if (tmp > 1)
      innerPartsSize.push_back(tmp);
     else
      break;
    }
   }
   partition2Level[++lClusterCnt] = originalHeight + 1;
   innerPartsSize.push_back(1);//The last partition has one element

  } else {
   /*if(minLevelDist<=0)
    minLevelDist=2;//default parameter*/
   while (innerPartsTmp > 1) {
    while (accuSlackGroups[originalHeight - levelCut - 1] <= innerPartsTmp
           && levelCut < levelNo) {
     levelCut++;
     if (originalHeight - levelCut - 1 < 0 || originalHeight - levelCut - 1 >= levelNo)
      break;
    }
    //Ensures a certain number of level in each partition
    innerPartsSize.push_back(innerPartsTmp);
    int tmp = originalHeight - levelCut - minLevelDist;
    if (tmp > partition2Level[lClusterCnt] && tmp < levelNo) {
     //Due to tuning parameter we need this check
     partition2Level[++lClusterCnt] = tmp;
    }
    innerPartsTmp /= divRate;
    preLevelCut = levelCut;
    levelCut = 0; // starting from the root
   }
   partition2Level[++lClusterCnt] = originalHeight + 1;
   innerPartsSize.push_back(1);//The last partition has one element
   delete[]accuSlackGroups;
  }
#if 0
  for (int i1 = 0; i1 <= lClusterCnt; ++i1) {
   std::cout << partition2Level[i1] << ",";
  }
  for (int i1 = lClusterCnt + 1; i1 < 7; ++i1) {
   std::cout << partition2Level[i1] << ",";
  }
 //std::cout<<"\n";
#endif
  return lClusterCnt;
 }

 void makeSlackedLevelSet(int n, int clusterCnt, int *partition2Level, int originalHeight,
                          std::vector<std::vector<int>> &slackGroups, std::vector<std::vector<int>> &slackedLevelSet,
                          int *node2Level) {
  std::vector<int> removedIDX;
  //taking a node from the node closer to the root
  for (int k = clusterCnt; k > 0; --k) {
   for (int i = partition2Level[k - 1]; i < partition2Level[k]; ++i) {
    //for (int i = partition2Level[k]-1; i > partition2Level[k-1]; --i) {
    //for each node in the level of the kth leveled partition
    for (int j = 0; j < slackGroups[i].size(); ++j) {
     //Assign the slack node to the lowest level possible
     for (int l = 0; l < k - 1; ++l) {
      int curSlackedNode = slackGroups[i][j];
      assert(curSlackedNode < n && curSlackedNode >= 0);
      //slackGroups[i].erase(slackGroups[i].begin()+j);
      //make sure the slacked node can be shifted to p-level l
      int targetedL = partition2Level[l + 1] - 1;
      assert(l + 1 < originalHeight);
      if (node2Level[curSlackedNode] < targetedL) {
       slackedLevelSet[targetedL].push_back(curSlackedNode);
       //slackGroups[i].erase(slackGroups[i].begin()) ;
       //Add the index to remove it later
       removedIDX.push_back(j);
      }
     }
    }
    //Removing the marked nodes
    for (int m = 0; m < removedIDX.size(); --m) {
     int tmpIDX = removedIDX[m];
     for (int m0 = 0; m0 < m; ++m0) {
      tmpIDX--;
     }
     slackGroups[i].erase(slackGroups[i].begin() + tmpIDX);
    }
    removedIDX.erase(removedIDX.begin(), removedIDX.end());
   }
  }
  //assign the remaining slacked nodes where they already are in slackgroup.
  for (int i = 0; i < clusterCnt; ++i) {
   for (int k = partition2Level[i]; k < partition2Level[i + 1]; ++k) {
    for (int j = 0; j < slackGroups[k].size(); ++j) {
     slackedLevelSet[k].push_back(slackGroups[k][j]);
    }
    slackGroups[k].erase(slackGroups[k].begin(), slackGroups[k].end());
   }
  }

#if 1//def VERIFY
  for (int m = 0; m < originalHeight; ++m) {
   assert(slackGroups[m].size() == 0);
  }
#endif

 }
}