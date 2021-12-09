//
// Created by kazem on 4/1/19.
//
// This part of the code is obtained from OSQP repository and is only for Ruitz Scaling
#ifndef PROJECT_OSQP_UTILS_H
#define PROJECT_OSQP_UTILS_H

#include "nasoq/common/def.h"

namespace nasoq {

# define MIN_SCALING (1e-04) ///< Minimum scaling value
# define MAX_SCALING (1e+04) ///< Maximum scaling value

 void vec_mult_scalar(double *a, double sc, int n);

 double vec_mean(const double *a, int n);


 void mat_inf_norm_cols_sym_triu(const CSC *M, double *E);

 void mat_premult_diag(CSC *A, const double *d);

 void mat_postmult_diag(CSC *A, const double *d);

 void vec_ew_prod(const double *a, const double *b, double *c, int n);

 void vec_ew_recipr(const double *a, double *b, int n);


 void vec_ew_sqrt(double *a, int n);


 void vec_set_scalar(double *a, double sc, int n);

// Set values lower than threshold SCALING_REG to 1
 void limit_scaling(double *D, int n);

 void mat_inf_norm_cols(const CSC *M, double *E);

 void vec_ew_max_vec(const double *a, const double *b, double *c, int n);


 void mat_inf_norm_rows(const CSC *M, double *E);


 double vec_norm_inf(const double *v, int l);

/* multiply scalar to matrix */
 void mat_mult_scalar(CSC *A, double sc);


/**
 * Compute infinite norm of the colums of the KKT matrix without forming it
 *
 * The norm is stored in the vector v = (D, E)
 *
 * @param P        Cost matrix
 * @param A        Contraints matrix
 * @param D        Norm of columns related to variables
 * @param D_temp_A Temporary vector for norm of columns of A
 * @param E        Norm of columns related to constraints
 * @param n        Dimension of KKT matrix
 */
 void compute_inf_norm_cols_KKT(const CSC *P, const CSC *A,
                                CSC *B,
                                double *D, double *D_temp_A,
                                double *E, double *F, int n);

}
#endif //PROJECT_OSQP_UTILS_H
