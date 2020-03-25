//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_DFS_H
#define CHOLOPENMP_DFS_H

#include <vector>
#include "def.h"

/* non-recursive version for actual use */

int dfsC		/* return the new value of k */
  (
    int p,		/* start the DFS at a root node p */
    int k,		/* start the node numbering at k */
    int *Post,	/* Post ordering, modified on output */
    int *Head ,	/* Head [p] = youngest child of p; EMPTY on output */
    int *Next ,	/* Next [j] = sibling of j; unmodified */
    int *Pstack 	/* workspace of size n, undefined on input or output */
  )
{
 int j, phead ;

 /* put the root node on the stack */
 Pstack [0] = p ;
 phead = 0 ;

 /* while the stack is not empty, do: */
 while (phead >= 0)
 {
  /* grab the node p from top of the stack and get its youngest child j */
  p = Pstack [phead] ;
  j = Head [p] ;
  if (j == EMPTY)
  {
   /* all children of p ordered.  remove p from stack and order it */
   phead-- ;
   Post [k++] = p ;	/* order node p as the kth node */
  }
  else
  {
   /* leave p on the stack.  Start a DFS at child node j by putting
    * j on the stack and removing j from the list of children of p. */
   Head [p] = Next [j] ;
   Pstack [++phead] = j ;
  }
 }
 return (k) ;	/* the next node will be numbered k */
}

//From CSPARSE
/* depth-first-search of the graph of a matrix, starting at node j */
int dfs (int j, int *Gp, int *Gi, int top, int *xi, int *pstack, const int *pinv)
{
 int i, p, p2, done, jnew, head = 0 ;
 //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
 //Gp = G->p ; Gi = G->i ;
 xi [0] = j ;                /* initialize the recursion stack */

 while (head >= 0)
 {
  j = xi [head] ;         /* get j from the top of the recursion stack */
  jnew = pinv ? (pinv [j]) : j ;
  if (!CS_MARKED (Gp, j))
  {
   CS_MARK (Gp, j) ;       /* mark node j as visited */
   pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ;
  }
  done = 1 ;                  /* node j done if no unvisited neighbors */
  p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;
  for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
  {
   i = Gi [p] ;            /* consider neighbor node i */
   if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
   pstack [head] = p ;     /* pause depth-first search of node j */
   xi [++head] = i ;       /* start dfs at node i */
   done = 0 ;              /* node j is not done */
   break ;                 /* break, to start dfs (i) */
  }
  if (done)               /* depth-first search at node j is done */
  {
   head-- ;            /* remove j from the recursion stack */
   xi [--top] = j ;    /* and place in the output stack */
  }
 }
 return (top) ;
}


/* depth-first-search of the graph of a matrix, starting at node j
 * depth-first-search of the graph of a matrix, starting at node j
 * modified to report the intersections
 * marked[i]=0; not visited
 * marked[i]=1; visited now
 * marked[i]=-1; visited in previous CC
 * clashedNodes has the list of nodes clashed with visited nodes
 * in previous CCs*/
int dfs_CSC_CC (size_t  n, int j, int *Gp, int *Gi,
                int *marked, int top, int *xi,
                int *pstack,
                std::vector<int> &clashedNodes,
                const int *pinv) {
 int i, p, p2, done, jnew, head = 0 ;
 //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
 //Gp = G->p ; Gi = G->i ;
 xi [0] = j ;                /* initialize the recursion stack */

 while (head >= 0) {
  j = xi [head] ;         /* get j from the top of the recursion stack */
  jnew = pinv ? (pinv [j]) : j ;
  if (!marked[jnew])  {
   marked[jnew] = 1;       /* mark node j as visited */
   pstack [head] =  Gp [jnew] ;
  }
  if(marked[jnew]==-1){//visited in previous CCs
   clashedNodes.push_back(jnew);
  }
  done = 1 ;                  /* node j done if no unvisited neighbors */
  p2 = Gp [jnew+1] ;
  for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
  {
   i = Gi [p] ;            /* consider neighbor node i */
   if (marked[i]==-1)
    clashedNodes.push_back(i);//Another node from prev CCs.
   if (marked[i]) continue ;   /* skip visited node i */
   pstack [head] = p ;     /* pause depth-first search of node j */
   xi [++head] = i ;       /* start dfs at node i */
   done = 0 ;              /* node j is not done */
   break ;                 /* break, to start dfs (i) */
  }
  if (done)               /* depth-first search at node j is done */
  {
   head-- ;            /* remove j from the recursion stack */
   xi [--top] = j ;    /* and place in the output stack */
  }
 }
 return (top) ;
}


/* depth-first-search of the graph of a matrix, starting at node j
 * modified to report the intersections
 * marked[i]=0; not visited
 * marked[i]=1; visited now
 * marked[i]=-1; visited in previous CC*/
int dfs_BCSC_CC (size_t  n, int j, size_t *Gp, size_t *Gi_ptr, int *Gi,
              const int *blk2Col,const int *col2Blk, int *marked,
              int top, int *xi,
              int *pstack,
                 std::vector<int> &clashedNodes,
              const int *pinv){
 int i, p, done, jnew, head = 0 ;
 //bool *marked = new bool[n]();
 //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
 //Gp = G->p ; Gi = G->i ;
 xi [0] = blk2Col[j] ;                /* initialize the recursion stack */
 while (head >= 0){
  j = xi [head];         /* get j from the top of the recursion stack */
  jnew=col2Blk[j];
  if (!marked[jnew]){//not visited before at all
   marked[jnew] = 1;       /* mark node j as visited */
   pstack [head] = Gi_ptr[j] ;
  }
  if(marked[jnew]==-1){//visited in previous CCs
   clashedNodes.push_back(jnew);
  }
  done = 1 ;                  /* node j done if no unvisited neighbors */

  int supWdth = blk2Col[jnew+1] - blk2Col[jnew];
  int nxtCol = blk2Col[jnew+1];
  int st=pstack [head];
  int p2 = Gi_ptr[nxtCol] ;
  for (int p1 = st ; p1 < p2 ; p1++)  /* examine all neighbors of j */
  {
   i = Gi [p1];            /* consider neighbor node i */
   int tmp = col2Blk[i];
   supWdth = blk2Col[tmp+1] - blk2Col[tmp];
   if (marked[tmp]==-1)
    clashedNodes.push_back(tmp);//Another node from prev CCs.
   if (marked[tmp]) continue ;   /* skip visited node i */
   pstack [head] = p1;     /* pause depth-first search of node j */
   xi [++head] = i ;       /* start dfs at node i */
   done = 0 ;              /* node j is not done */
   break ;                 /* break, to start dfs (i) */
  }
  if (done)               /* depth-first search at node j is done */
  {
   head-- ;            /* remove j from the recursion stack */
   xi [--top] = col2Blk[j] ;    /* and place in the output stack */
  }
 }
 return (top) ;
}


//FIXME needs testing
/* depth-first-search of the graph of a matrix, starting at node j */
int dfs_BCSC (size_t  n, int j, size_t *Gp, size_t *Gi_ptr, int *Gi,
              const int *blk2Col,const int *col2Blk, bool *marked,
              int top, int *xi,
              int *pstack,
              const int *pinv){
 int i, p, done, jnew, head = 0 ;
 //bool *marked = new bool[n]();
 //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
 //Gp = G->p ; Gi = G->i ;
 xi [0] = blk2Col[j] ;                /* initialize the recursion stack */
 while (head >= 0){
  j = xi [head];         /* get j from the top of the recursion stack */
  jnew=col2Blk[j];
  if (!marked[jnew]){
   marked[jnew] = true;       /* mark node j as visited */
   pstack [head] = Gi_ptr[j] ;
  }
  done = 1 ;                  /* node j done if no unvisited neighbors */

  int supWdth = blk2Col[jnew+1] - blk2Col[jnew];
  int st=pstack [head];
  int p2 = Gi_ptr[j+supWdth] ;
  for (int p1 = st ; p1 < p2 ; p1++)  /* examine all neighbors of j */
  {
   i = Gi [p1];            /* consider neighbor node i */
   int tmp = col2Blk[i];
   supWdth = blk2Col[tmp+1] - blk2Col[tmp];
   if (marked[tmp]) continue ;   /* skip visited node i */
   pstack [head] = p1;     /* pause depth-first search of node j */
   xi [++head] = i ;       /* start dfs at node i */
   done = 0 ;              /* node j is not done */
   break ;                 /* break, to start dfs (i) */
  }
  if (done)               /* depth-first search at node j is done */
  {
   head-- ;            /* remove j from the recursion stack */
   xi [--top] = col2Blk[j] ;    /* and place in the output stack */
  }
 }
 return (top) ;
}

/* depth-first search and postorder of a tree rooted at node j */
int tdfs(int j, int k, int *head, const int *next, int *post, int *stack)
{
 int i, p, top = 0 ;
 if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
 stack [0] = j ;                 /* place j on the stack */
 while (top >= 0)                /* while (stack is not empty) */
 {
  p = stack [top] ;           /* p = top of stack */
  i = head [p] ;              /* i = youngest child of p */
  if (i == -1)
  {
   top-- ;                 /* p has no unordered children left */
   post [k++] = p ;        /* node p is the kth postordered node */
  }
  else
  {
   head [p] = next [i] ;   /* remove i from children of p */
   stack [++top] = i ;     /* start dfs on child node i */
  }
 }
 return (k) ;
}


#endif //CHOLOPENMP_DFS_H
