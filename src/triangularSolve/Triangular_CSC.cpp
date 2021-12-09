//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/triangularSolve/Triangular_CSC.h"

#include <chrono>

#include "nasoq/common/Reach.h"

namespace nasoq {

 int lsolve(int n, int *Lp, int *Li, double *Lx, double *x) {
  int p, j;
  if (!Lp || !Li || !x) return (0);                     /* check inputs */
  for (j = 0; j < n; j++) {
   x[j] /= Lx[Lp[j]];
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
    x[Li[p]] -= Lx[p] * x[j];
   }
  }
  return (1);
 }

 int ltsolve(int n, int *Lp, int *Li, double *Lx, double *x) {
  int p, j;
  if (!Lp || !Li || !x) return (0);                      /* check inputs */
  for (j = n - 1; j >= 0; j--) {
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
    x[j] -= Lx[p] * x[Li[p]];
   }
   x[j] /= Lx[Lp[j]];
  }
  return (1);
 }

 unsigned long flopCoutLSolve(int n, int *Lp, int *Li, double *Lx, double *x) {
  int p, j;
  unsigned long flopCount = 0;
  if (!Lp || !Li || !x) return (0);                     /* check inputs */
  for (j = 0; j < n; j++) {
   x[j] /= Lx[Lp[j]];
   flopCount++;
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
    x[Li[p]] -= Lx[p] * x[j];
    flopCount += 2;
   }
  }
  return (flopCount);
 }

 int lsolvePar(int n, int *Lp, int *Li, double *Lx, double *x, int levels, int *levelPtr, int *levelSet, int chunk) {
  if (!Lp || !Li || !x) return (0);                     /* check inputs */
  for (int l = 0; l < levels; ++l) {
   int li = 0;
#pragma omp parallel for \
   default(shared) private(li)  \
   schedule(static)
   for (li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
    int j = levelSet[li];
    x[j] /= Lx[Lp[j]];
    for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
     double tmp = Lx[p] * x[j];
     int idx = Li[p];
#pragma omp atomic
     x[idx] -= tmp;
    }
   }
  }
  return (1);
 }

 int lsolveParH2(int n, int *Lp, int *Li, double *Lx, double *x, int levels, int *levelPtr, int *levelSet, int parts,
                 int *parPtr, int *partition, int chunk) {
  if (!Lp || !Li || !x) return (0);                     /* check inputs */
  for (int i1 = 0; i1 < levels; ++i1) {
#pragma omp parallel //shared(lValues)//private(map, contribs)
   {
#pragma omp  for schedule(static)
    for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
     for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
      int j = partition[k1];
      x[j] /= Lx[Lp[j]];
      //    #pragma omp critical
      for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
       double tmp = Lx[p] * x[j];
       int idx = Li[p];
#pragma omp atomic
       x[idx] -= tmp;
      }
     }
    }
   }
  }
  return (1);
 }

 int lsolvePar2(int n, int *Lp, int *Li, double *Lx, double *x) {
  int p, j;
  if (!Lp || !Li || !x) return (0);                     /* check inputs */
  for (j = 0; j < n; j++) {
   x[j] /= Lx[Lp[j]];
   //#pragma omp parallel for
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
    x[Li[p]] -= Lx[p] * x[j];
   }
  }
  return (1);
 }

 int lsolve_reach_dec(int n, int *Gp, int *Gi, double *Gx, int *Bp, int *Bi, double *Bx, int k, int *xi, double *x,
                      const int *pinv, double &symDuration) {
  int j, p, px, top;
  std::chrono::time_point<std::chrono::system_clock> start, end;
  top = reach(n, Gp, Gi, Bp, Bi, k, xi, pinv);
  start = std::chrono::system_clock::now();
  for (px = top; px < n; px++) {
   j = xi[px];
   x[j] /= Gx[(Gp[j])];
   p = Gp[j] + 1;
   for (; p < Gp[j + 1]; p++) {
    x[Gi[p]] -= Gx[p] * x[j];
   }
  }
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> tmp = end - start;
  symDuration = tmp.count();
  return (top);
 }

 void lSolveSympiler(int n, int *Lp, int *Li, const double *Lx, double *x, int *reachSet, int reachSetSize) {
  int p, px, j;
  x[0] /= Lx[0]; // Peel col 0
  double x_Li_1 = Lx[1] * x[0];
  double x_Li_2 = Lx[2] * x[0];
  x[*(Li + 1)] -= x_Li_1;
  x[*(Li + 2)] -= x_Li_2;

  for (px = 1; px < 3; px++) {
   j = reachSet[px];
   x[j] /= Lx[Lp[j]];
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
    x[Li[p]] -= Lx[p] * x[j];
  }
  x[7] /= Lx[20]; // Peel col 7
  double x_Li_21 = Lx[21] * x[7];
  double x_Li_22 = Lx[22] * x[7];
  x[*(Li + 21)] -= x_Li_21;
  x[*(Li + 22)] -= x_Li_22;

  for (px = 4; px < reachSetSize; px++) {
   j = reachSet[px];
   x[j] /= Lx[Lp[j]];
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
    x[Li[p]] -= Lx[p] * x[j];
  }
 }
}
