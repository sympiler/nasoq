//
// Created by kazem on 4/21/19.
//

#ifndef PROJECT_UPDOWN_TEST_H
#define PROJECT_UPDOWN_TEST_H

#include <vector>

#include "nasoq/common/def.h"

namespace nasoq {

/*
* Build super matrix from A and B,  All are stored in CSC.
 * C is eq
 * B is ineq
*/
 int build_super_matrix_with_eq(CSC *A, CSC *B, CSC *C, double reg_diag,
                                size_t &SM_size, size_t &SM_nz,
                                int *&SMp, int *&SMi, double *&SMx,
                                double *&sm_rhs);


/*
 * A is inequality matrix
 */
 void build_super_solve_with_eq(CSC *H, const CSC *A, CSC *B, const double *rhs,
                                double reg_diag, const std::vector<int> &mod_col, double *x_all,
                                int outer_iter = 2, int inner_iter = 2, double stop_tol = 1e-15);


 void build_super_solve_with_eq_mkl(CSC *H, const CSC *A, CSC *B,
                                    const double *rhs, double reg_diag,
                                    const std::vector<int> &mod_col, double *x_all,
                                    int outer_iter = 2, int inner_iter = 2, double stop_tol = 1e-15);
}
#endif //PROJECT_UPDOWN_TEST_H
