//
// Created by kazem on 3/11/19.
//

#ifndef PROJECT_QP_UTILS_H
#define PROJECT_QP_UTILS_H

#include <chrono>
#include <string>
#include <vector>

#include "nasoq/common/def.h"

namespace nasoq {
 double dot1(int n, double *a, double *b);


 int build_qp_from_file(std::string hessian_file,
                        std::string ineq_file,
                        std::string ineq_file_rhs,
                        std::string linear_file,
                        std::string opt_primal_file,
                        std::string opt_dual_file,
                        std::string opt_obj_file,
                        std::string warm_start_file,
                        size_t &sizeH, size_t &nnzH, double *&q,
                        int *&colH, int *&rowH, double *&valH,
                        size_t &ad1, size_t &ad2,
                        size_t &nnzA,
                        int *&colA, int *&rowA, double *&valA,
                        double *&b_ineq, double *&optimal_primal,
                        double &objective, double *&optimal_dual,
                        double *&init_x);

 int build_qp_optimality_from_file(std::string hessian_file,
                                   std::string ineq_file,
                                   std::string ineq_file_rhs,
                                   std::string linear_file,
                                   std::string opt_primal_file,
                                   std::string opt_dual_file,
                                   size_t &sizeH, size_t &nnzH, double *&q,
                                   int *&colH, int *&rowH, double *&valH,
                                   size_t &ad1, size_t &ad2,
                                   size_t &nnzA,
                                   int *&colA, int *&rowA, double *&valA,
                                   double *&b_ineq, double *&optimal_primal,
                                   double *&optimal_dual);

/*
 * takes 1 iterations, just add
 */
 void build_qp_01(size_t &sizeH, size_t &nnzH, double *&q,
                  int *&colH, int *&rowH, double *&valH,
                  size_t &ad1, size_t &ad2,
                  size_t &nnzA,
                  int *&colA, int *&rowA, double *&valA,
                  double *&b_ineq, double *&optimal_primal,
                  double &objective);

/*
 * takes 1 iterations, just add
 */
 void build_qp_02(size_t &sizeH, size_t &nnzH, double *&q,
                  int *&colH, int *&rowH, double *&valH,
                  size_t &ad1, size_t &ad2,
                  size_t &nnzA,
                  int *&colA, int *&rowA, double *&valA,
                  double *&b_ineq, double *&optimal_primal,
                  double &objective);

/*
 * takes 2 iterations, just add
 */
 void build_qp_03(size_t &sizeH, size_t &nnzH, double *&q,
                  int *&colH, int *&rowH, double *&valH,
                  size_t &ad1, size_t &ad2,
                  size_t &nnzA,
                  int *&colA, int *&rowA, double *&valA,
                  double *&b_ineq, double *&optimal_primal,
                  double &objective);

/*
 * eq constraints
 */
 void build_qp_04(size_t &sizeH, size_t &nnzH, double *&q,
                  int *&colH, int *&rowH, double *&valH,
                  size_t &ad1, size_t &ad2,
                  size_t &nnzA,
                  int *&colA, int *&rowA, double *&valA,
                  double *&b_ineq,
                  size_t &bd1, size_t &bd2,
                  size_t &nnzB,
                  int *&colB, int *&rowB, double *&valB,
                  double *&b_eq,
                  double *&optimal_primal,
                  double &objective);

 struct qp_info {
  std::chrono::time_point<std::chrono::system_clock> tot_st, tot_end;
  std::chrono::time_point<std::chrono::system_clock> fct_st, fct_end;
  std::chrono::time_point<std::chrono::system_clock> ini_st, ini_end;
  bool sw;
  std::chrono::duration<double> elapsed_seconds;
  double tot;
  double factt, init;
  int num_update, num_downdate, num_refactor, num_solve;

  qp_info();

  std::chrono::time_point<std::chrono::system_clock> tic();

  std::chrono::time_point<std::chrono::system_clock> toc();

  double elapsed_time(std::chrono::time_point<std::chrono::system_clock> beg,
                      std::chrono::time_point<std::chrono::system_clock> lst);


  void print();
 };

 struct idx_val {
  double val;
  int idx;

  idx_val();

  idx_val(int i, double v);
 };

 bool CMP(idx_val a, idx_val b);

 int remove_from_list(std::vector<idx_val> &list_const, int c_no);


/*
 * Computes primal objective
 */
 double compute_primal_obj(double *primal_vars, CSC *H, double *q);


/*
* Computes min(b-Bx,0)
*/
 double constraint_sat_norm(CSC *B, CSC *A, double *b, double *a,
                            double *primal_vars);


/*
* Computes the norm of Hx + q + B^T*dual_vars
*/
 double lagrangian_residual_norm(CSC *H, CSC *B, CSC *BT, CSC *A, CSC *AT,
                                 double *q, double *primal_vars,
                                 double *dual_vars, double *dual_vars_eq);


/*
 * Computes ||dual_vars * b-Bx||inf
 * ||dual_vars_eq * a-Ax||inf
 */
 double complementarity_norm(CSC *B, CSC *A, double *b, double *a,
                             double *primal_vars, double *dual_vars);

 double non_negativity_norm(CSC *B, double *dual_vars);

/*
 * The Fischer-Burmeister (FB) function is f(a,b) = a + b - sqrt(a^2 + b^2)
 */
 void Fischer_Burmeister_func(int n, double *a, double *b, double *f);

/*
 * Converting bound constraint l < x < u to Cx < d
 */
 void bound_to_general_ineq();
}
#endif //PROJECT_QP_UTILS_H
