//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/symbolic/Partitioning.h"

namespace nasoq {

 int partitioning(int inSize, int *inTree, int *inCost, int *inChildPtr, int *inChildNo, int *nChild, int n,
                  int partitionNum, int &outSize, int *outCost, int *outNode2Par,
                  std::vector<std::vector<int>> &parList) {
  auto *visited = new bool[inSize]();
  auto *mark = new bool[inSize]();
  int Threshold, k = 0;
  Threshold = n / partitionNum;//Almost the same numjber of columns in each partition
  //Threshold=1;
  std::vector<int> stack;
  //Initial partitioning, each partition one node.
  //parList.resize(inSize);
  std::vector<int> extraDim;
  parList.push_back(extraDim);
#if 0
  for (int l = 0; l < inSize; ++l) {
         std::cout<<l<<": "<<inTree[l]<<",";
     }
     std::cout<<"\n";
#endif
  int curPart = 0;
  outSize = 0;
  for (int curNode = 0; curNode < inSize; ++curNode) {//k is in partitioned node
   if (inTree[curNode] == -2)
    continue;
   if (!visited[curNode]) {
    stack.push_back(curNode);
    while (stack.size() != 0) {
     k = stack[0];
     if (nChild[k] == 0) {
      //Add it to current part
      stack.erase(stack.begin());
      parList[curPart].push_back(k);
      outCost[curPart] += inCost[k];
      visited[k] = true;//mark as visited
      //nChild[inTree[k]] =k>=0? nChild[inTree[k]]-1:nChild[inTree[k]];
      if (inTree[k] >= 0) {//The node is a single node
       nChild[inTree[k]]--;
       if (!visited[inTree[k]]) {
        if (!mark[inTree[k]])//if it is not in the stack
         stack.insert(stack.begin(), inTree[k]);//Add its parent
        mark[inTree[k]] = true;//There is at least a children of this node in this partition
       }
      } else//There is no other node to check, check next one
       break;
     } else {
      //adding its children for further evaluation
//                    std::cout<<"The children of: "<<k<<"\n";
      for (int i = inChildPtr[k]; i < inChildPtr[k + 1]; ++i) {
       int tmpChild = inChildNo[i];
       if (!visited[tmpChild]) {
        if (nChild[tmpChild] == 0) {//First add leaves
         //    stack.insert(stack.begin(),tmpChild);
         //    tmpFront++;
         if (!visited[k]) {
          int tmpk = k;
          if (stack.size() > 0) {
           while (tmpk != stack[stack.size() - 1]) {//Mark all parents upwards, they should be in the partiotion
            mark[tmpk] = true;//There is at least a children of this node in this partition
            tmpk = inTree[tmpk];
           }
          }
          mark[tmpk] = true;
         }
         parList[curPart].push_back(tmpChild);
         outCost[curPart] += inCost[tmpChild];
         visited[tmpChild] = true;//mark as visited
         nChild[inTree[tmpChild]]--;//we know inTree[tmpChild] is k
        } else {
         stack.insert(stack.begin(), tmpChild);
        }
	   }
      }
     }//End else
     if (outCost[curPart] > Threshold) {
      for (int s = 0; s < stack.size(); ++s) {
       int leftover = stack[s];
       if (mark[leftover]) {
        parList[curPart].push_back(leftover);
        outCost[curPart] += inCost[leftover];
        visited[leftover] = true;//mark as visited
        mark[leftover] = false;//for further iterations
        if (inTree[leftover] >= 0)
         nChild[inTree[leftover]]--;
       }
      }
      stack.erase(stack.begin(), stack.end());
      break;
     }
    }
    //It is either reaches the defined size or there is no other node to add (stack is empty)
    curPart++;
    std::vector<int> extraDim;
    parList.push_back(extraDim);
   }
  }
  //remove the last empty element
  parList.erase(parList.begin() + parList.size() - 1);
  outSize = curPart;
  for (int i = 0; i < curPart; ++i) {
   for (int j = 0; j < parList[i].size(); ++j) {
    outNode2Par[parList[i][j]] = i;
   }
  }

#if 0
  for (int i = 0; i < curPart; ++i) {
         std::cout<<"Partition #"<<i<<" :";
         for (int j = 0; j < parList[i].size(); ++j) {
             std::cout<<parList[i][j]<<",";
         }
         std::cout<<"\n";
     }
#endif
#if 0
  auto *check = new bool[inSize]();
     for (int i = 0; i < curPart; ++i) {
         for (int j = 0; j < parList[i].size(); ++j) {
                 check[parList[i][j]]=true;
         }
     }
     for (int i = 0; i < inSize; ++i) {
         assert(check[i]);
     }
     delete []check;
#endif
  delete[]mark;
  delete[]visited;
  return outSize;
 }
}