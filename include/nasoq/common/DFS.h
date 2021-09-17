//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_DFS_H
#define CHOLOPENMP_DFS_H

#include <vector>
#include <cstdlib>
namespace nasoq {
/* non-recursive version for actual use */

 int dfsC  /* return the new value of k */
   (
     int p,  /* start the DFS at a root node p */
     int k,  /* start the node numbering at k */
     int *Post, /* Post ordering, modified on output */
     int *Head, /* Head [p] = youngest child of p; EMPTY on output */
     int *Next, /* Next [j] = sibling of j; unmodified */
     int *Pstack  /* workspace of size n, undefined on input or output */
   );

//From CSPARSE
/* depth-first-search of the graph of a matrix, starting at node j */
 int dfs(int j, int *Gp, int *Gi, int top, int *xi, int *pstack, const int *pinv);


/* depth-first-search of the graph of a matrix, starting at node j
 * depth-first-search of the graph of a matrix, starting at node j
 * modified to report the intersections
 * marked[i]=0; not visited
 * marked[i]=1; visited now
 * marked[i]=-1; visited in previous CC
 * clashedNodes has the list of nodes clashed with visited nodes
 * in previous CCs*/
 int dfs_CSC_CC(size_t n, int j, int *Gp, int *Gi,
                int *marked, int top, int *xi,
                int *pstack,
                std::vector<int> &clashedNodes,
                const int *pinv);


/* depth-first-search of the graph of a matrix, starting at node j
 * modified to report the intersections
 * marked[i]=0; not visited
 * marked[i]=1; visited now
 * marked[i]=-1; visited in previous CC*/
 int dfs_BCSC_CC(size_t n, int j, size_t *Gp, size_t *Gi_ptr, int *Gi,
                 const int *blk2Col, const int *col2Blk, int *marked,
                 int top, int *xi,
                 int *pstack,
                 std::vector<int> &clashedNodes,
                 const int *pinv);


//FIXME needs testing
/* depth-first-search of the graph of a matrix, starting at node j */
 int dfs_BCSC(size_t n, int j, size_t *Gp, size_t *Gi_ptr, int *Gi,
              const int *blk2Col, const int *col2Blk, bool *marked,
              int top, int *xi,
              int *pstack,
              const int *pinv);

/* depth-first search and postorder of a tree rooted at node j */
 int tdfs(int j, int k, int *head, const int *next, int *post, int *stack);

}
#endif //CHOLOPENMP_DFS_H
