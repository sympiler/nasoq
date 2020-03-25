//
// Created by kazem on 7/8/18.
//

#ifndef PROJECT_BFS_H
#define PROJECT_BFS_H

#include <vector>
#include <cstddef>
#include <assert.h>
/*
 * Mixed BFS and topological sort for BCSC
 */
int modifiedBFS(int n,
                size_t *Lp,
                size_t *Li_ptr,
                int* Li,
                const int* sup2col,
                const int* col2sup,
                int *inDegree,
                bool *visited,
                int *node2partition,
                int* &levelPtr,
                size_t* levelSet,
                int bfsLevel,
                std::vector<std::vector<int>> &newLeveledParList){

  std::vector<int> queue;
 //Let's do BFS for every leaf node of the
 for (int ii = levelPtr[bfsLevel]; ii < levelPtr[bfsLevel+1]; ++ii) {
  int curNode=levelSet[ii];
  assert(node2partition[curNode]>=0);
  queue.push_back(curNode);
  while (!queue.empty()){
   int popedNode=queue[0];
   queue.erase(queue.begin());
   visited[popedNode]=true;
   newLeveledParList[node2partition[popedNode]].push_back(popedNode);
   //Find the adjacent nodes
   int curCol = sup2col[popedNode];
   int nxtCol = sup2col[popedNode+1];
   int supWdt = nxtCol-curCol;
   for (int r = Li_ptr[curCol]+supWdt; r < Li_ptr[nxtCol]; ++r) {
    int cn=col2sup[Li[r]];
    inDegree[cn]--;
    if(inDegree[cn]==1 && !visited[cn]){
     queue.push_back(cn);
    }
   }
  }
 }
}


/*
 * Mixed BFS and topological sort for CSC format.
 */
int modifiedBFS_CSC(int n,
                int *Lp,
                int* Li,
                int *inDegree,
                bool *visited,
                int *node2partition,
                int* &levelPtr,
                int* levelSet,
                int bfsLevel,
                std::vector<std::vector<int>> &newLeveledParList){

 std::vector<int> queue;
 //Let's do BFS for every leaf node of the
 for (int ii = levelPtr[bfsLevel]; ii < levelPtr[bfsLevel+1]; ++ii) {
  int curNode=levelSet[ii];
  assert(node2partition[curNode]>=0);
  queue.push_back(curNode);
  while (!queue.empty()){
   int popedNode=queue[0];
   queue.erase(queue.begin());
   visited[popedNode]=true;
   newLeveledParList[node2partition[popedNode]].push_back(popedNode);
   //Find the adjacent nodes
   for (int r = Lp[popedNode]; r < Lp[popedNode+1]; ++r) {
    int cn=Li[r];
    inDegree[cn]--;
    if(inDegree[cn]==1 && !visited[cn]){
     queue.push_back(cn);
    }
   }
  }
 }

}
#endif //PROJECT_BFS_H
