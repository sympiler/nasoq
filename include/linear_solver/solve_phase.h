//
// Created by kazem on 11/23/18.
//

#ifndef PROJECT_SOLVE_PHASE_H
#define PROJECT_SOLVE_PHASE_H

#include <cstddef>

namespace nasoq{
/*
 * Scale a vector by another vector vec2 = vec2/vec1
 */
void scale_vec_vec(int n, double *vec1, double *vec2);

/*
 * Cholesky non-supernodal solve phase
 */
void solve_phase_simplicial_ll(int n, int* Lp, int* Li, double* Lx, double *x);

/*
 * LDL non-supernodal solve phase
 */
void solve_phase_simplicial_ldl(int n, int* Lp, int* Li, double* Lx,
                               double *d_val,
                               double *x);

/*
 * Solve phase for LLT factorization, static pivoting
 * the input is a blocked L (BCSC)
 * rhs is in input and the solution will be stored in the output
 */
void solve_phase_ll_blocked(size_t n, double *x, int *col2sup, int *sup2col,
                             size_t *newCol, int *newRow, double *newVal,
                             size_t *rowP, size_t nBlocks,
                             int nnz);

/*
 * Solve phase for LDLT factorization, static pivoting
 * the input is a blocked L (BCSC)
 */
void solve_phase_ldl(size_t n, double *d_val, double *x, int *col2sup,
                     int *sup2col,
                     size_t *newCol, int *newRow, double *newVal,
                     size_t *rowP, size_t nBlocks,
                     int nnz);

/*
 * Solve phase for LDLT factorization, static pivoting
 * the input is a blocked L (BCSC)
 */
void solve_phase_ldl_blocked(size_t n, double *d_val, double *x, int *col2sup, int *sup2col,
                     size_t *newCol, int *newRow, double *newVal, size_t *rowP, size_t nBlocks,
                     int nnz);


/*
 * Solve phase for LDLT factorization
 * the input is a blocked L (BCSC)
 */
void solve_phase_ldl_blocked_parallel(size_t n, double *d_val, double *x, int *col2sup, int *sup2col,
                                      size_t *newCol, int *newRow, double *newVal,
                                      size_t *rowP, size_t nBlocks,int nnz,
                                      int levels, int *levelPtr, int *levelSet,
                                      int parts,  int *parPtr, int *partition,
                                      int chunk
                             );


/*
 * Solve phase for LDLT factorization
 * the input is a blocked L (BCSC)
 */
void solve_phase_ldl_blocked_parallel_permuted(size_t n, double *d_val, double *x,
                                      int *col2sup, int *sup2col,
                                      size_t *newCol, int *newRow, double *newVal,
                                      size_t *rowP, size_t nBlocks,int nnz,
                                      int levels, int *levelPtr, int *levelSet,
                                      int parts,  int *parPtr, int *partition, int chunk,
                                      int *perm, int *iperm);


/*
 * Solve phase for LDLT factorization
 * the input is a blocked L (BCSC)
 */
void solve_phase_ldl_blocked_parallel_permuted_update
(size_t n, double *d_val, double *x, int *col2sup, int *sup2col,
  size_t *newCol, int *newRow, double *newVal,
  size_t *rowP, size_t nBlocks,int nnz,
  int levels, int *levelPtr, int *levelSet,
  int parts,  int *parPtr, int *partition, int chunk,
  int s_level_no, int *s_level_ptr, int *s_level_set,
  bool *mask,
  int *mask_col, double *ws=NULL);



}
#endif //PROJECT_SOLVE_PHASE_H
