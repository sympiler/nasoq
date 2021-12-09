//
// Created by kazem on 3/16/19.
//

#ifndef PROJECT_SOLVER_GENERATION_H
#define PROJECT_SOLVER_GENERATION_H

#include "nasoq/common/Util.h"

namespace nasoq {
 int build_linear_solve_from_file(std::string matrix_in,
                                  std::string rhs_file,
                                  std::string opt_file,
                                  size_t &sizeH, size_t &nnzH, double *&rhs,
                                  int *&colH, int *&rowH, double *&valH,
                                  double *&solution) {

  CSC *H_tmp = new CSC;
  if (!readMatrix(matrix_in, H_tmp->ncol, H_tmp->nzmax, H_tmp->p,
                  H_tmp->i, H_tmp->x))
   return 0;
  //print_csc("dd", H_tmp->ncol, H_tmp->p, H_tmp->i, H_tmp->x);
  bool is_expanded = expandMatrix(H_tmp->ncol, H_tmp->nzmax, H_tmp->p,
                                  H_tmp->i, H_tmp->x,
                                  nnzH, colH,
                                  rowH, valH);
  if (is_expanded) {
   //H has the expanded version already.
   sizeH = H_tmp->ncol;
   allocateAC(H_tmp, 0, 0, 0, FALSE);
  } else {
   valH = H_tmp->x;
   colH = H_tmp->p;
   rowH = H_tmp->i;
   nnzH = H_tmp->nzmax;
   sizeH = H_tmp->ncol;
  }
  nnzH = colH[sizeH];


  if (rhs_file != "none") {
   rhs = new double[sizeH];
   read_vector(rhs_file, sizeH, rhs);
   if (opt_file != "none") {
    solution = new double[sizeH];
    read_vector(opt_file, sizeH, solution);
   } else {
    solution = NULL;
   }
  } else {//Generate RHS such that solution==1
   rhs = new double[sizeH];
   solution = new double[sizeH];
   rhsInit_linearSolver_onlylower(sizeH, colH, rowH, valH, rhs);
   for (int i = 0; i < sizeH; ++i) {
    solution[i] = 1.0;
   }
  }
  return 1;
 }


 int build_linear_solve_update_from_file(std::string matrix_in,
                                         std::string rhs_file,
                                         std::string update_in,
                                         std::string update_rhs_in,
                                         std::string opt_file,
                                         size_t &sizeH, size_t &nnzH,
                                         double *&rhs,
                                         int *&colH, int *&rowH, double *&valH,
                                         size_t &row_A, size_t &col_A, size_t &nnzA,
                                         double *&a, int *&colA, int *&rowA,
                                         double *&valA,
                                         double *&solution) {

  if (!readMatrix(matrix_in, sizeH, nnzH, colH, rowH, valH))
   return 0;
  if (!readMatrix_rect(update_in, row_A, col_A, nnzA, colA, rowA, valA))
   return 0;
  a = new double[row_A];
  if (!read_vector(update_rhs_in, row_A, a))
   return 0;

  solution = NULL;
  if (rhs_file != "none") {
   rhs = new double[sizeH];
   read_vector(rhs_file, sizeH, rhs);
   if (opt_file != "none") {
    solution = new double[sizeH];
    read_vector(opt_file, sizeH, solution);
   }
  } else {//Generate RHS such that solution==1
   rhs = new double[sizeH];
   solution = new double[sizeH];
   rhsInit_linearSolver_onlylower(sizeH, colH, rowH, valH, rhs);
   for (int i = 0; i < sizeH; ++i) {
    solution[i] = 1.0;
   }
  }
  return 1;
 }
}
#endif //PROJECT_SOLVER_GENERATION_H
