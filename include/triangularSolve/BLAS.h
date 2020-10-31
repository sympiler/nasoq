//
// Created by kazem on 7/18/17.
//

#ifndef TRIANGOPENMP_BLAS_H
#define TRIANGOPENMP_BLAS_H

namespace nasoq {

 void dlsolve_blas_nonUnit(int ldm, int ncol, double *M, double *rhs);

 void lSolve_dense_col_sync(int colSize, int col, double *M, double *rhs);

 void dmatvec_blas(
   int ldm, /* in -- leading dimension of M */
   int nrow, /* in */
   int ncol, /* in */
   double *M, /* in */
   double *vec, /* in */
   double *Mxvec /* in/out */
 );

}
#endif //TRIANGOPENMP_BLAS_H
