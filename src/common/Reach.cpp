//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/common/Reach.h"

#include <chrono>

#include "nasoq/common/def.h"
#include "nasoq/common/DFS.h"

namespace nasoq {

 int reach(int n, int *Gp, int *Gi, int *Bp, int *Bi, int k, int *xi, const int *pinv) {
  int p, top;
  if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1);    /* check inputs */
  //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
  top = n;
  for (p = Bp[k]; p < Bp[k + 1]; p++) {
   if (!CS_MARKED (Gp, Bi[p]))    /* start a dfs at unmarked node i */
   {
    top = dfs(Bi[p], Gp, Gi, top, xi, xi + n, pinv);
   }
  }
  for (p = top; p < n; p++) CS_MARK (Gp, xi[p]);  /* restore G */
  return (top);
 }

 int reach_sn(int n, int *Gp, int *Gi, int *Bp, int *Bi, int k, int *PBset, const int *pinv, int sn, int *col2sup,
              double &analysis) {
  int p, top, PBsize = 0;
  bool *checked = new bool[sn]();
  int *xi = new int[2 * n]();
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  start = std::chrono::system_clock::now();
  if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1);    /* check inputs */
  //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
  top = n;
  for (p = Bp[k]; p < Bp[k + 1]; p++) {
   if (!CS_MARKED (Gp, Bi[p]))    /* start a dfs at unmarked node i */
   {
    top = dfs(Bi[p], Gp, Gi, top, xi, xi + n, pinv);
   }
  }
  for (p = top; p < n; p++) CS_MARK (Gp, xi[p]);  /* restore G */
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  analysis = elapsed_seconds.count();
  for (int i = top; i < n; ++i) {
   if (!checked[col2sup[xi[i]]]) {
    checked[col2sup[xi[i]]] = true;
    PBset[PBsize++] = col2sup[xi[i]];
   }
  }
#if DEBUG > 0
  for(int i=0; i<PBsize; i++){
         std::cout<<PBset[i]<<",";
     }
     std::cout<<"\n";
#endif
  delete[]checked;
  delete[]xi;
  return (PBsize);
 }

 int reach_col(int n, int *Gp, int *Gi, int *Bp, int *Bi, int k, int *xi, const int *pinv) {
  int p, top;
  if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1);    /* check inputs */
  //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
  top = n;
  for (p = Bp[k]; p < Bp[k + 1]; p++) {
   if (!CS_MARKED (Gp, Bi[p]))    /* start a dfs at unmarked node i */
   {
    top = dfs(Bi[p], Gp, Gi, top, xi, xi + n, pinv);
   }
  }
  for (p = top; p < n; p++) CS_MARK (Gp, xi[p]);  /* restore G */
  return (top);
 }

 int ereach(int n, int *Ap, int *Ai, int k, const int *parent, int *s, int *w) {
  int i, p, len, top;
  if (!Ap || !Ai || !parent || !s || !w) return (-1);   /* check inputs */
  top = n;
  CS_MARK (w, k);                /* mark node k as visited */
  for (p = Ap[k]; p < Ap[k + 1]; p++) {
   i = Ai[p];                /* A(i,k) is nonzero */
   if (i > k) continue;       /* only use upper triangular part of A */
   for (len = 0; !CS_MARKED (w, i); i = parent[i]) /* traverse up etree*/
   {
    s[len++] = i;         /* L(k,i) is nonzero */
    CS_MARK (w, i);        /* mark i as visited */
   }
   while (len > 0) s[--top] = s[--len]; /* push path onto stack */
  }
  for (p = top; p < n; p++) CS_MARK (w, s[p]);    /* unmark all nodes */
  CS_MARK (w, k);                /* unmark node k */
  return (top);                  /* s [top..n-1] contains pattern of L(k,:)*/
 }

 int ereach_sn(int n, int *Ap, int *Ai, int col1, int col2, int *col2sup, const int *parent, int *s, int *w) {
  int i, p, len, top;
  if (!Ap || !Ai || !parent || !s || !w) return (-1);   /* check inputs */
  top = n;
  for (int k = col1; k < col2; ++k) {
   ASSERT(col2sup[k] < n);
   if (k == col1) CS_MARK (w, col2sup[k]);                /* mark node k as visited */
   for (p = Ap[k]; p < Ap[k + 1]; p++) {
    i = col2sup[Ai[p]];                /* A(i,k) is nonzero block */
    ASSERT(i < n);
    if (Ai[p] > k)
     continue;       /* only use upper triangular part of A */
    //if(col2sup[i] == col2sup[Ai[p-1]]) continue; // from the same supernode
    for (len = 0; !CS_MARKED (w, i) && i >= 0 && parent[i] != INVISIBLE;
         i = parent[i]) /* traverse up etree*/
    {
     s[len++] = i;         /* L(k,i) is nonzero */
     ASSERT(len < n);
     CS_MARK (w, i);        /* mark i as visited */
     ASSERT(i < n);
    }
    while (len > 0) s[--top] = s[--len]; /* push path onto stack */
   }
  }
  for (p = top; p < n; p++) CS_MARK (w, s[p]);    /* unmark all nodes */
  //for (int k = col1; k < col2; ++k) {
  CS_MARK (w, col2sup[col1]);                /* unmark node k */
  //}
  return (top);                  /* s [top..n-1] contains pattern of L(k,:)*/
 }
}