//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/triangularSolve/BLAS.h"

namespace nasoq {

 void dlsolve_blas_nonUnit(int ldm, int ncol, double *M, double *rhs)//general triangular solver
 {
  int k;
  double x0, x1, x2, x3, x4, x5, x6, x7;
  double *M0;
  double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
  int firstcol = 0;

  M0 = &M[0];

  while (firstcol < ncol - 7) { /* Do 8 columns */
   Mki0 = M0;
   Mki1 = Mki0 + ldm + 1;
   Mki2 = Mki1 + ldm + 1;
   Mki3 = Mki2 + ldm + 1;
   Mki4 = Mki3 + ldm + 1;
   Mki5 = Mki4 + ldm + 1;
   Mki6 = Mki5 + ldm + 1;
   Mki7 = Mki6 + ldm + 1;

   x0 = rhs[firstcol] / *Mki0++;
   x1 = (rhs[firstcol + 1] - x0 * *Mki0++) / *Mki1++;
   x2 = (rhs[firstcol + 2] - x0 * *Mki0++ - x1 * *Mki1++) / *Mki2++;
   x3 = (rhs[firstcol + 3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++) / *Mki3++;
   x4 = (rhs[firstcol + 4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
         - x3 * *Mki3++) / *Mki4++;
   x5 = (rhs[firstcol + 5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
         - x3 * *Mki3++ - x4 * *Mki4++) / *Mki5++;
   x6 = (rhs[firstcol + 6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
         - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++) / *Mki6++;
   x7 = (rhs[firstcol + 7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
         - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
         - x6 * *Mki6++) / *Mki7++;

   rhs[firstcol++] = x0;
   rhs[firstcol++] = x1;
   rhs[firstcol++] = x2;
   rhs[firstcol++] = x3;
   rhs[firstcol++] = x4;
   rhs[firstcol++] = x5;
   rhs[firstcol++] = x6;
   rhs[firstcol++] = x7;

   for (k = firstcol; k < ncol; k++)
    rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
             - x2 * *Mki2++ - x3 * *Mki3++
             - x4 * *Mki4++ - x5 * *Mki5++
             - x6 * *Mki6++ - x7 * *Mki7++;

   M0 += 8 * ldm + 8;
  }

  while (firstcol < ncol - 3) { /* Do 4 columns */
   Mki0 = M0;
   Mki1 = Mki0 + ldm + 1;
   Mki2 = Mki1 + ldm + 1;
   Mki3 = Mki2 + ldm + 1;

   x0 = rhs[firstcol] / *Mki0++;
   x1 = (rhs[firstcol + 1] - x0 * *Mki0++) / *Mki1++;
   x2 = (rhs[firstcol + 2] - x0 * *Mki0++ - x1 * *Mki1++) / *Mki2++;
   x3 = (rhs[firstcol + 3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++) / *Mki3++;

   rhs[firstcol++] = x0;
   rhs[firstcol++] = x1;
   rhs[firstcol++] = x2;
   rhs[firstcol++] = x3;

   for (k = firstcol; k < ncol; k++)
    rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
             - x2 * *Mki2++ - x3 * *Mki3++;

   M0 += 4 * ldm + 4;
  }

  if (firstcol < ncol - 1) { /* Do 2 columns */
   Mki0 = M0;
   Mki1 = Mki0 + ldm + 1;

   x0 = rhs[firstcol] / *Mki0++;
   x1 = (rhs[firstcol + 1] - x0 * *Mki0++) / *Mki1++;

   rhs[firstcol++] = x0;
   rhs[firstcol++] = x1;

   for (k = firstcol; k < ncol; k++)
    rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
   M0 += 2 * ldm + 2;
  }

  if (firstcol == ncol - 1) { /* Do 1 columns */
   Mki0 = M0;
   x0 = rhs[firstcol] / *Mki0;
   rhs[firstcol] = x0;
  }
 }

 void lSolve_dense_col_sync(int colSize, int col, double *M, double *rhs) {
//#pragma omp critical
  for (int i = 0; i < col; ++i) {
//#pragma omp atomic
   rhs[i] /= M[i * colSize + i];
   for (int j = i + 1; j < col; ++j) {
    double tmp = M[i * colSize + j] * rhs[i];
//#pragma omp atomic
    rhs[j] -= tmp;
   }
  }
  //return 1;
 }

 void dmatvec_blas(int ldm, int nrow, int ncol, double *M, double *vec, double *Mxvec) {
  double vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
  double *M0;
  double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
  int firstcol = 0;
  int k;

  M0 = &M[0];
  while (firstcol < ncol - 7) { /* Do 8 columns */

   Mki0 = M0;
   Mki1 = Mki0 + ldm;
   Mki2 = Mki1 + ldm;
   Mki3 = Mki2 + ldm;
   Mki4 = Mki3 + ldm;
   Mki5 = Mki4 + ldm;
   Mki6 = Mki5 + ldm;
   Mki7 = Mki6 + ldm;

   vi0 = vec[firstcol++];
   vi1 = vec[firstcol++];
   vi2 = vec[firstcol++];
   vi3 = vec[firstcol++];
   vi4 = vec[firstcol++];
   vi5 = vec[firstcol++];
   vi6 = vec[firstcol++];
   vi7 = vec[firstcol++];

   for (k = 0; k < nrow; k++)
    Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
                + vi2 * *Mki2++ + vi3 * *Mki3++
                + vi4 * *Mki4++ + vi5 * *Mki5++
                + vi6 * *Mki6++ + vi7 * *Mki7++;

   M0 += 8 * ldm;
  }

  while (firstcol < ncol - 3) { /* Do 4 columns */

   Mki0 = M0;
   Mki1 = Mki0 + ldm;
   Mki2 = Mki1 + ldm;
   Mki3 = Mki2 + ldm;

   vi0 = vec[firstcol++];
   vi1 = vec[firstcol++];
   vi2 = vec[firstcol++];
   vi3 = vec[firstcol++];
   for (k = 0; k < nrow; k++)
    Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
                + vi2 * *Mki2++ + vi3 * *Mki3++;

   M0 += 4 * ldm;
  }

  while (firstcol < ncol) {  /* Do 1 column */

   Mki0 = M0;
   vi0 = vec[firstcol++];
   for (k = 0; k < nrow; k++)
    Mxvec[k] += vi0 * *Mki0++;

   M0 += ldm;
  }

 }
}