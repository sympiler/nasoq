//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/PostOrder.h"

#include "nasoq/common/def.h"
#include "nasoq/common/DFS.h"
#include "nasoq/common/SparseUtils.h"

namespace nasoq {

 int postOrderC(int *Parent, size_t n, int *Weight, int *Post, int status) {
  int *Head, *Next, *Pstack, *Iwork;
  int j, p, k, w, nextj;
  size_t s;
  int ok = TRUE;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  /* RETURN_IF_NULL_COMMON (EMPTY) ;
   RETURN_IF_NULL (Parent, EMPTY) ;
   RETURN_IF_NULL (Post, EMPTY) ;
   Common->status = CHOLMOD_OK ;*/

  /* ---------------------------------------------------------------------- */
  /* allocate workspace */
  /* ---------------------------------------------------------------------- */

  /* s = 2*n */
  s = mult_size_t(n, 3, &ok);
  if (!ok) {
//        ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
   return (EMPTY);
  }

  //CHOLMOD(allocate_work) (n, s, 0, Common) ;
  /*if (Common->status < CHOLMOD_OK)
  {
      return (EMPTY) ;
  }*/
  //ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */
  Iwork = new int[s + 1]();
//    Head  = Common->Head ;	/* size n+1, initially all EMPTY */
//    Iwork = Common->Iwork ;
  Head = Iwork;
  Next = Iwork + n + 1;  /* size n (i/i/l) */
  Pstack = Iwork + 2 * n + 1; /* size n (i/i/l) */

  for (int i = 0; i < n + 1; ++i) {
   Head[i] = EMPTY;
  }
  /* ---------------------------------------------------------------------- */
  /* construct a link list of children for each node */
  /* ---------------------------------------------------------------------- */

  if (Weight == NULL) {

   /* in reverse order so children are in ascending order in each list */
   for (j = n - 1; j >= 0; j--) {
    p = Parent[j];
    if (p >= 0 && p < ((int) n)) {
     /* add j to the list of children for node p */
     Next[j] = Head[p];
     Head[p] = j;
    }
   }

   /* Head [p] = j if j is the youngest (least-numbered) child of p */
   /* Next [j1] = j2 if j2 is the next-oldest sibling of j1 */

  } else {

   /* First, construct a set of link lists according to Weight.
    *
    * Whead [w] = j if node j is the first node in bucket w.
    * Next [j1] = j2 if node j2 follows j1 in a link list.
    */

   int *Whead = Pstack;     /* use Pstack as workspace for Whead [ */

   for (w = 0; w < ((int) n); w++) {
    Whead[w] = EMPTY;
   }
   /* do in forward order, so nodes that ties are ordered by node index */
   for (j = 0; j < ((int) n); j++) {
    p = Parent[j];
    if (p >= 0 && p < ((int) n)) {
     w = Weight[j];
     w = MAX (0, w);
     w = MIN (w, ((int) n) - 1);
     /* place node j at the head of link list for weight w */
     Next[j] = Whead[w];
     Whead[w] = j;
    }
   }

   /* traverse weight buckets, placing each node in its parent's list */
   for (w = n - 1; w >= 0; w--) {
    for (j = Whead[w]; j != EMPTY; j = nextj) {
     nextj = Next[j];
     /* put node j in the link list of its parent */
     p = Parent[j];
     ASSERT (p >= 0 && p < ((int) n));
     Next[j] = Head[p];
     Head[p] = j;
    }
   }

   /* Whead no longer needed ] */
   /* Head [p] = j if j is the lightest child of p */
   /* Next [j1] = j2 if j2 is the next-heaviest sibling of j1 */
  }

  /* ---------------------------------------------------------------------- */
  /* start a DFS at each root node of the etree */
  /* ---------------------------------------------------------------------- */

  k = 0;
  for (j = 0; j < ((int) n); j++) {
   if (Parent[j] == EMPTY) {
    /* j is the root of a tree; start a DFS here */
    k = dfsC(j, k, Post, Head, Next, Pstack);
   }
  }

  /* this would normally be EMPTY already, unless Parent is invalid */
  for (j = 0; j < ((int) n); j++) {
   Head[j] = EMPTY;
  }

/*    PRINT1 (("postordered "ID" nodes\n", k)) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;*/
  delete[]Iwork;
  return (k);
 }

 int *postOrder(const int *parent, int n) {
/* post order a forest
 * Obtained from CSparse library
 * */
  int j, k = 0, *post, *w, *head, *next, *stack;
  if (!parent) return (NULL);                        /* check inputs */
  //post = cs_malloc (n, sizeof (csi)) ;                /* allocate result */
  //w = cs_malloc (3*n, sizeof (csi)) ;                 /* get workspace */
  post = new int[n];
  w = new int[3 * n];
  //if (!w || !post) return (cs_idone (post, NULL, w, 0)) ;
  if (!w || !post)
   return NULL;
  head = w;
  next = w + n;
  stack = w + 2 * n;
  for (j = 0; j < n; j++) head[j] = -1;           /* empty linked lists */
  for (j = n - 1; j >= 0; j--)            /* traverse nodes in reverse order*/
  {
   if (parent[j] == -1) continue;    /* j is a root */
   next[j] = head[parent[j]];      /* add j to list of its parent */
   head[parent[j]] = j;
  }
  for (j = 0; j < n; j++) {
   if (parent[j] != -1) continue;    /* skip j if it is not a root */
   k = tdfs(j, k, head, next, post, stack);
  }
  // return (cs_idone (post, NULL, w, 1)) ;  /* success; free w, return post */
  delete[]w;
  return post;
 }
}