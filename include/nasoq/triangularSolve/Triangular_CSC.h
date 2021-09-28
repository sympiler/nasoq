//
// Created by kazem on 7/18/17.
//

#ifndef TRIANGOPENMP_TRIANGULAR_CSC_H
#define TRIANGOPENMP_TRIANGULAR_CSC_H

namespace nasoq {
/*
 ****** Serial implementation
 */
 int lsolve(int n, int *Lp, int *Li, double *Lx, double *x);

/*
 * L^T x = b
 */
 int ltsolve(int n, int *Lp, int *Li, double *Lx, double *x);

/*
 * Counting the number of FLOPS in triangular solve
 */
 unsigned long flopCoutLSolve(int n, int *Lp, int *Li, double *Lx, double *x);


/*
 ****** Parallel
 */
 int lsolvePar(int n, int *Lp, int *Li, double *Lx, double *x,
               int levels, int *levelPtr, int *levelSet, int chunk);

/*
 ****** Parallel H2
 */
 int lsolveParH2(int n, int *Lp, int *Li, double *Lx, double *x,
                 int levels, int *levelPtr, int *levelSet,
                 int parts, int *parPtr, int *partition,
                 int chunk);


/*
 *
 */
 int lsolvePar2(int n, int *Lp, int *Li, double *Lx, double *x);


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
 int lsolve_reach_dec(int n, int *Gp, int *Gi, double *Gx, int *Bp, int *Bi, double *Bx,
                      int k, int *xi, double *x, const int *pinv, double &symDuration);

/*
 * only for motive example
 */
 void lSolveSympiler(int n, int *Lp, int *Li,
                     const double *Lx, double *x,
                     int *reachSet, int reachSetSize);
}
#endif //TRIANGOPENMP_TRIANGULAR_CSC_H
