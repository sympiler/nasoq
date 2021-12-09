//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/linear_solver/solve_phase.h"

#include "nasoq/common/Sym_BLAS.h"
#include "nasoq/triangularSolve/Triangular_BCSC.h"
#include "nasoq/triangularSolve/Triangular_CSC.h"

namespace nasoq {

 void scale_vec_vec(int n, double *vec1, double *vec2) {
  for (int i = 0; i < n; ++i) {
   double tmp = *vec2 / *(vec1++);
   *(vec2++) = tmp;
   //*(vec2++) = *vec2 / *(vec1++);
  }
 }

 void solve_phase_simplicial_ll(int n, int *Lp, int *Li, double *Lx, double *x) {
  lsolve(n,Lp,Li,Lx,x);
  //print_vec("\nlsolve: ",0,n,x);
  ltsolve(n,Lp,Li,Lx,x);
  //print_vec("\nltsolve: ",0,n,x);
 }

 void solve_phase_simplicial_ldl(int n, int *Lp, int *Li, double *Lx, double *d_val, double *x) {
  lsolve(n,Lp,Li,Lx,x);
  //print_vec("\nlsolve: ",0,n,x);
  //***** diagonal scaling
  scale_vec_vec(n, d_val, x);
  ltsolve(n,Lp,Li,Lx,x);
  //print_vec("\nltsolve: ",0,n,x);
 }

 void
 solve_phase_ll_blocked(size_t n, double *x, int *col2sup, int *sup2col, size_t *newCol, int *newRow, double *newVal,
                        size_t *rowP, size_t nBlocks, int nnz) {
  //*************** Serial Blocked
  blockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
  //*************** Serial Blocked FWD solve
  blockedLTsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
 }

 void solve_phase_ldl(size_t n, double *d_val, double *x, int *col2sup, int *sup2col, size_t *newCol, int *newRow,
                      double *newVal, size_t *rowP, size_t nBlocks, int nnz) {//*************** Serial Blocked
  blockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
#if 0
  for (int i = 0; i < n; ++i) {
  std::cout<<i<<"->"<<x[i]<<"/"<<d_val[i]<<"= "<<x[i]/d_val[i]<<"\n";
 }
#endif
  //***** diagonal scaling
  scale_vec_vec(n, d_val, x);
#if 0
  for (int i = 0; i < n; ++i) {
  std::cout<<i<<"->"<<x[i]<<"\n";
 }
#endif
  //*************** Serial Blocked FWD solve
  blockedLTsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);

 }

 void
 solve_phase_ldl_blocked(size_t n, double *d_val, double *x, int *col2sup, int *sup2col, size_t *newCol, int *newRow,
                         double *newVal, size_t *rowP, size_t nBlocks, int nnz) {//*************** Serial Blocked
  blockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
#if 0
  for (int i = 0; i < n; ++i) {
  std::cout<<i<<"->"<<x[i]<<"/"<<d_val[i]<<"= "<<x[i]/d_val[i]<<"\n";
 }
#endif
  //***** diagonal scaling
  //scale_vec_vec(n, d_val, x);
  blocked_2by2_solver(n,d_val,x,1,1,n);
#if 0
  for (int i = 0; i < n; ++i) {
  std::cout<<i<<"->"<<x[i]<<"\n";
 }
#endif
  //*************** Serial Blocked FWD solve
  blockedLTsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);

 }

 void solve_phase_ldl_blocked_parallel(size_t n, double *d_val, double *x, int *col2sup, int *sup2col, size_t *newCol,
                                       int *newRow, double *newVal, size_t *rowP, size_t nBlocks, int nnz, int levels,
                                       int *levelPtr, int *levelSet, int parts, int *parPtr, int *partition, int chunk) {
  //*************** Parallel L solve
  H2LeveledBlockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                         sup2col, nBlocks, x, levels, levelPtr, levelSet,
                         parts, parPtr, partition, chunk);
#if 0
  for (int i = 0; i < n; ++i) {
  std::cout<<i<<"->"<<x[i]<<"/"<<d_val[i]<<"= "<<x[i]/d_val[i]<<"\n";
 }
#endif
  //***** diagonal scaling
  blocked_2by2_solver(n,d_val,x,1,1,n);
#if 0
  for (int i = 0; i < n; ++i) {
  std::cout<<i<<"->"<<x[i]<<"\n";
 }
#endif
  //*************** Parallel Blocked BWD solve
  H2LeveledBlockedLTsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                          sup2col, nBlocks, x, levels, levelPtr, levelSet,
                          parts, parPtr, partition, chunk);
 }

 void solve_phase_ldl_blocked_parallel_permuted(size_t n, double *d_val, double *x, int *col2sup, int *sup2col,
                                                size_t *newCol, int *newRow, double *newVal, size_t *rowP,
                                                size_t nBlocks, int nnz, int levels, int *levelPtr, int *levelSet,
                                                int parts, int *parPtr, int *partition, int chunk, int *perm,
                                                int *iperm) {
  //*************** Parallel L solve
  H2LeveledBlockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                         sup2col, nBlocks, x, levels, levelPtr, levelSet,
                         parts, parPtr, partition, chunk);
  //***** diagonal scaling
  blocked_2by2_solver(n,d_val,x,1,1,n);
//*************** Parallel Blocked BWD solve
  H2LeveledBlockedLTsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                          sup2col, nBlocks, x, levels, levelPtr, levelSet,
                          parts, parPtr, partition, chunk);
 }

 void solve_phase_ldl_blocked_parallel_permuted_update(size_t n, double *d_val, double *x, int *col2sup, int *sup2col,
                                                       size_t *newCol, int *newRow, double *newVal, size_t *rowP,
                                                       size_t nBlocks, int nnz, int levels, int *levelPtr,
                                                       int *levelSet, int parts, int *parPtr, int *partition, int chunk,
                                                       int s_level_no, int *s_level_ptr, int *s_level_set, bool *mask,
                                                       int *mask_col, double *ws) {
  if(s_level_no<0){
   //*************** Parallel L solve
   H2LeveledBlockedLsolve_update(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                                 sup2col, nBlocks, x, levels, levelPtr, levelSet,
                                 parts, parPtr, partition, chunk,mask,ws);
   //***** diagonal scaling
   blocked_2by2_solver_update(n,d_val,x,1,1,n,mask_col);
   //*************** Parallel Blocked BWD solve
   H2LeveledBlockedLTsolve_update(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                                  sup2col, nBlocks, x, levels, levelPtr, levelSet,
                                  parts, parPtr, partition, chunk,mask,ws);
  }else{
   leveledBlockedLsolve_update(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                               sup2col, nBlocks, x,
                               s_level_no, s_level_ptr, s_level_set,
                               chunk, mask);
   //***** diagonal scaling
   blocked_2by2_solver_update(n,d_val,x,1,1,n,mask_col);
//*************** Parallel Blocked BWD solve
   LeveledBlockedLTsolve_update(n, newCol, newRow, newVal, nnz, rowP, col2sup,
                                sup2col, nBlocks, x, s_level_no, s_level_ptr,
                                s_level_set,
                                chunk, mask,ws);
  }

 }
}