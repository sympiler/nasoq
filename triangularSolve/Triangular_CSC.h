//
// Created by kazem on 7/18/17.
//

#ifndef TRIANGOPENMP_TRIANGULAR_CSC_H
#define TRIANGOPENMP_TRIANGULAR_CSC_H

#include <immintrin.h>
#include "../common/Reach.h"

/*
 ****** Serial implementation
 */
int lsolve (int n, int* Lp, int* Li, double* Lx, double *x){
 int p, j;
 if (!Lp || !Li || !x) return (0) ;                     /* check inputs */
 for (j = 0 ; j < n ; j++)
 {
  x [j] /= Lx [Lp [j]] ;
  for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
  {
   x [Li [p]] -= Lx [p] * x [j] ;
  }
 }
 return (1) ;
}

/*
 * L^T x = b
 */
int ltsolve (int n, int* Lp, int* Li, double* Lx, double *x)
{
 int p, j;
 if (!Lp || !Li || !x) return (0) ;                      /* check inputs */
 for (j = n-1 ; j >= 0 ; j--)
 {
  for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
  {
   x [j] -= Lx [p] * x [Li [p]] ;
  }
  x [j] /= Lx [Lp [j]] ;
 }
 return (1) ;
}

/*
 * Counting the number of FLOPS in triangular solve
 */
unsigned long flopCoutLSolve (int n, int* Lp, int* Li, double* Lx, double *x){
 int p, j;
 unsigned long flopCount=0;
 if (!Lp || !Li || !x) return (0) ;                     /* check inputs */
 for (j = 0 ; j < n ; j++){
  x [j] /= Lx [Lp [j]] ;
  flopCount++;
  for (p = Lp [j]+1 ; p < Lp [j+1] ; p++){
   x [Li [p]] -= Lx [p] * x [j] ;
   flopCount+=2;
  }
 }
 return (flopCount) ;
}


/*
 ****** Parallel
 */
int lsolvePar (int n, int* Lp, int* Li, double* Lx, double *x,
               int levels, int *levelPtr, int *levelSet, int chunk){
 if (!Lp || !Li || !x) return (0) ;                     /* check inputs */
 for (int l = 0; l < levels; ++l) {
  int li=0;
#pragma omp parallel for \
   default(shared) private(li)  \
   schedule(auto)
  for ( li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
   int j = levelSet[li];
   x [j] /= Lx [Lp [j]] ;
   for (int p = Lp [j]+1 ; p < Lp [j+1] ; p++)
   {
    double tmp =  Lx [p] * x [j] ;
    int idx = Li[p];
#pragma omp atomic
    x [idx] -= tmp ;
   }
  }
 }
 return (1) ;
}

/*
 ****** Parallel H2
 */
int lsolveParH2 (int n, int* Lp, int* Li, double* Lx, double *x,
               int levels, int *levelPtr, int *levelSet,
                 int parts,  int *parPtr, int *partition,
                 int chunk){
 if (!Lp || !Li || !x) return (0) ;                     /* check inputs */
 for (int i1 = 0; i1 < levels ; ++i1) {
#pragma omp parallel //shared(lValues)//private(map, contribs)
  {
#pragma omp  for schedule(auto)
   for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
    for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
     int j = partition[k1];
     x[j] /= Lx[Lp[j]];
 //    #pragma omp critical
     for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
      double tmp = Lx[p] * x[j];
      int idx = Li[p];
#pragma omp atomic
      x[idx] -= tmp ;
     }
    }
   }
  }
 }
 return (1) ;
}


/*
 *
 */
int lsolvePar2 (int n, int* Lp, int* Li, double* Lx, double *x){
 int p, j;
 if (!Lp || !Li || !x) return (0) ;                     /* check inputs */
 for (j = 0 ; j < n ; j++)
 {
  x [j] /= Lx [Lp [j]] ;
  //#pragma omp parallel for
  for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
  {
   x [Li [p]] -= Lx [p] * x [j] ;
  }
 }
 return (1) ;
}


/*
 * Vectorized implementation
 */
#if 0
typedef union
{
 __m256d v;
 double d[4];
} v4df_t;

int lsolveVectorize(int n, int* Lp, int* Li, const double* Lx, double *x) {
 double xx;
 v4df_t reg_Lx;
 v4df_t reg_x;
 v4df_t result0, result1, result2, result3;
 int mod=0;

#if 0
 for (int k = st ; k < bd1inReach ; k++)
 {
  j = reach[k];
  xx=x [j];
  xx /= Lx [Lp [j]] ;
  for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
  {
   x [Li [p]] -= Lx [p] * xx;
  }
  x[j]=xx;
 }
#endif

#if 0
 for (int k = bd1 ; k < bd2 ; k++) {
        xx = x[k];
        xx /= Lx[Lp [k]];
        reg_x.v = _mm256_set1_pd(xx);
        mod = (Lp [k+1] - Lp [k] - 1) % 4;
        for (int i1 = Lp [k] + 1; i1 < Lp [k+1] - mod; i1 += 4) {
            reg_Lx.v = _mm256_load_pd((double *) (Lx + i1));
            result0.v = _mm256_mul_pd(reg_Lx.v, reg_x.v);

            x[Li[i1]] -= result0.d[0];
            x[Li[i1 + 1]] -= result0.d[1];
            x[Li[i1 + 2]] -= result0.d[2];
            x[Li[i1 + 3]] -= result0.d[3];
        }
        for (int i1 = Lp [k+1] - mod; i1 < Lp [k+1]; ++i1) {
            x[Li[i1]] -= Lx[i1] * xx;
        }
        x[k] = xx;
    }
#endif
#if 1
 for (int k = 0 ; k < n ; k++) {
  xx = x[k];
  xx /= Lx[Lp [k]];
  if(xx != 0){
   reg_x.v = _mm256_set1_pd(xx);
   mod = (Lp [k+1] - Lp [k] - 1) % 16;
   for (int i1 = Lp [k] + 1; i1 < Lp [k+1] - mod; i1 += 16) {
    reg_Lx.v = _mm256_load_pd((double *) (Lx + i1));
    result0.v = _mm256_mul_pd(reg_Lx.v, reg_x.v);
    reg_Lx.v = _mm256_load_pd((double *) (Lx + i1+4));
    result1.v = _mm256_mul_pd(reg_Lx.v, reg_x.v);
    reg_Lx.v = _mm256_load_pd((double *) (Lx + i1+8));
    result2.v = _mm256_mul_pd(reg_Lx.v, reg_x.v);
    reg_Lx.v = _mm256_load_pd((double *) (Lx + i1+12));
    result3.v = _mm256_mul_pd(reg_Lx.v, reg_x.v);

    x[Li[i1]] -= result0.d[0];
    x[Li[i1 + 1]] -= result0.d[1];
    x[Li[i1 + 2]] -= result0.d[2];
    x[Li[i1 + 3]] -= result0.d[3];
    x[Li[i1 + 4]] -= result1.d[0];
    x[Li[i1 + 5]] -= result1.d[1];
    x[Li[i1 + 6]] -= result1.d[2];
    x[Li[i1 + 7]] -= result1.d[3];
    x[Li[i1 + 8]] -= result2.d[0];
    x[Li[i1 + 9]] -= result2.d[1];
    x[Li[i1 + 10]] -= result2.d[2];
    x[Li[i1 + 11]] -= result2.d[3];
    x[Li[i1 + 12]] -= result3.d[0];
    x[Li[i1 + 13]] -= result3.d[1];
    x[Li[i1 + 14]] -= result3.d[2];
    x[Li[i1 + 15]] -= result3.d[3];
   }
   for (int i1 = Lp [k+1] - mod; i1 < Lp [k+1]; ++i1) {
    x[Li[i1]] -= Lx[i1] * xx;
   }
   x[k] = xx;

  }
 }
#endif

}
#endif


/*
 * Pruned
 */
int lsolve_reach_dec (int n, int* Gp, int* Gi, double* Gx, int* Bp, int* Bi, double* Bx,
                      int k, int *xi, double *x, const int *pinv, double &symDuration)
{
 int j, p, px, top ;
 std::chrono::time_point<std::chrono::system_clock> start, end;
 top = reach (n, Gp, Gi, Bp, Bi, k, xi, pinv) ;
 start = std::chrono::system_clock::now();
 for (px = top ; px < n ; px++){
  j = xi [px] ;
  x [j] /= Gx [(Gp [j])] ;
  p = Gp [j]+1 ;
  for ( ; p < Gp [j+1] ; p++)
  {
   x [Gi [p]] -= Gx [p] * x [j] ;
  }
 }
 end = std::chrono::system_clock::now();
 std::chrono::duration<double>tmp=end-start;
 symDuration=tmp.count();
 return (top) ;
}

/*
 * only for motive example
 */
int lSolveSympiler(int n, int* Lp, int* Li,
                    const double* Lx, double *x,
                   int *reachSet, int reachSetSize){
 int p,px,j;
 x[0] /= Lx[0]; // Peel col 0
 double x_Li_1 = Lx[1] * x[0];
 double x_Li_2 = Lx[2] * x[0];
 x[*(Li+1)] -= x_Li_1;
 x[*(Li+2)] -= x_Li_2;

 for(px=1;px<3;px++){
  j=reachSet[px];
  x[j]/=Lx[Lp[j]];
  for(p=Lp[j]+1;p<Lp[j+1];p++)
   x[Li[p]]-=Lx[p]*x[j];
 }
 x[7] /= Lx[20]; // Peel col 7
 double x_Li_21 = Lx[21] * x[7];
 double x_Li_22 = Lx[22] * x[7];
 x[*(Li+21)] -= x_Li_21;
 x[*(Li+22)] -= x_Li_22;

 for(px=4;px<reachSetSize;px++){
  j=reachSet[px];x[j]/=Lx[Lp[j]];
  for(p=Lp[j]+1;p<Lp[j+1];p++)
   x[Li[p]]-=Lx[p]*x[j];
 }
}

#endif //TRIANGOPENMP_TRIANGULAR_CSC_H
