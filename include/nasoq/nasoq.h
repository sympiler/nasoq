//
// Created by kazem on 2/26/19.
//

#ifndef PROJECT_QP_DUAL_H
#define PROJECT_QP_DUAL_H
//#define CHOLROWMOD

#include <string>
#include <vector>

#include "nasoq/common/def.h"
#include "nasoq/QP/qp_utils.h"
#include "nasoq/QP/linear_solver_wrapper.h"

namespace nasoq {
#ifdef CHOLROWMOD
#include "cholmod_utils.h"
 //using namespace nasoq_cholmod;
#endif

/*
 * Class for passing setting to a QP solver. It is a superset of all
 * important input parameters of a QP solver.
 */
 struct QPSettings {
  double eps; // one threshold for all
  double eps_primal, eps_dual, eps_slack, eps_nn; // If thresholds are requested per KKT conditions
  double diag_perturb;
  int batch_size;
  double eps_rel;
  int scaling;
  double zero_thresh;
  int inner_iter_ref;
  int outer_iter_ref;
  int max_iter;
  double stop_tol;
  int max_iter_nas;
  std::string nasoq_variant;

  QPSettings();
 };

 enum nasoq_mode {
  Fixed = 0, AUTO, Tuned, PREDET
 };
 enum nasoq_status {
  Optimal = 1, Inaccurate = 2, NotConverged = 3, Infeasible=0, NotFinished=4
 };

 struct nasoq_config {
  int inner_iter, outer_iter;
  double pert_diag, stop_tol;

  nasoq_config(int a, int b, double c, double d);
 };

/*
* min 0.5*xHx + qx s.t. Ax=a and Bx≤b
*/
 struct Nasoq {
  std::string sol_name;
  int hessian_size, eq_const_size, ineq_const_size;
  CSC *H, *HT, *A, *AT, *B, *BT;

  double *q, *a, *b;
  int n_active;
  double diag_perturb, zero_thresh, eps_abs, eps_rel;
  double primal_obj, dual_obj, objective;
  double non_negativity_infn, complementarity_infn;
  int inner_iter_ref, outer_iter_ref, max_iter;
  double stop_tol;

/// Symbolic info
  int *eTree; //skkt_size
  int *col2const, *const2col;
  int *used_const;

/// Super KKT
  CSC *sKKT, *sKKTt;
  double *sKKTrhs;
  int *etree;
  int *pinv;
  int *col2sup;
  size_t skkt_col;

  int scaling;
  double *D, *E, *F, *Dinv, *Einv, *Finv;
  double c, cinv;

  double *rec_length;

  double *dual_vars, *dual_vars_eq; //size of ineq_const_size
  double *primal_vars;
//double *x0; // Initial solution, warm-start
  double *kkt_solution, *lagrange_mult, *lagrange_mult_eq, *descent;
  std::vector<int> active_set;
  int *as0;// Initial active set, warm-start
  int is_converged; // 0/1
  int num_iter, max_iter_nas;
  int used_constraints, num_active; //
  int auto_reg_en;
  int warm_start;
  int batch_size;
  double *workspace;
  double qp_scalar, inv_qp_scalar;
  nasoq_mode variant;
  nasoq_status ret_val;
  SolverSettings *ss;

/// Profiling
  qp_info *qi;
  double *dual_FB_a, *dual_FB_b, *dual_FB;
  double *primal_FB_a, *primal_FB_b, *primal_FB;
  double lag_res, cons_sat_norm;

#ifdef CHOLROWMOD
  cholmod_factor *L ;
  cholmod_sparse *kkt;
  CSC *kkt_updated;
  cholmod_common Common, *cm;
  size_t k_size;
#endif

  Nasoq(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in,
        size_t B_row, size_t B_col, int *Bp, int *Bi,
        double *Bx, double *b_ineq);

  Nasoq(size_t H_size, int *Hp, int *Hi, double *Hx, double *q_in,
        size_t A_size1, size_t A_size2, int *Ap, int *Ai, double *Ax,
        double *a_eq,
        size_t B_size1, size_t B_size2, int *Bp, int *Bi, double *Bx,
        double *b_ineq);

  ~Nasoq();

/*
* Setting default setting for QP
*/
  void default_setting();

#ifdef CHOLROWMOD
  int symbolic_QP_cholmod(){
   int status =  0;
   transpose_unsym(B->nrow,B->ncol,B->p,B->i,B->x,
                   BT->nrow,BT->ncol,BT->p,BT->i,BT->x);
   BT->nzmax=B->nzmax;
   //print_csc("BT:\n",BT->ncol,BT->p,BT->i,BT->x);
   //build_super_kkt();
   //Creates a new solver instance
   if(A->nrow>0){
    transpose_unsym(A->nrow,A->ncol,A->p,A->i,A->x,
                    AT->nrow,AT->ncol,AT->p,AT->i,AT->x);
    AT->nzmax=A->nzmax;
   }

   /*ss->ldl_variant =4;
   ss->ldl_update_variant=2;
   ss->solver_mode = 1;
   ss->req_ref_iter=outer_iter_ref;
   if(outer_iter_ref>0 && inner_iter_ref ==0){
    ss->max_inner_iter = 1; // will be wrong if it is zero
   }else{
    ss->max_inner_iter = inner_iter_ref;
   }*/

   size_t total_nnz = H->nzmax + B->nzmax;
   double min_v,max_v;
   //max_min_spmat(H->ncol, H->p, H->i, H->x,max_v,min_v);
   //setting perturbation
   if(auto_reg_en==2){
    double H_norm2 = norm_sparse_int(H->ncol, H->p, H->i, H->x, -1, 2);
    auto_perturbation2(H->ncol, total_nnz,H_norm2, max_v-min_v,
                       B->nrow, B->nrow/double(H->ncol), eps_abs,
                       diag_perturb,
                       outer_iter_ref,inner_iter_ref,stop_tol);
   } else if(auto_reg_en==1){
    determine_iterations(total_nnz, eps_abs,outer_iter_ref,inner_iter_ref);
    stop_tol = 1e-15;
   }
    diag_perturb = pow(10, -6);
   /*ss->tol_abs = ss->tol_rel = stop_tol;
   ss->req_ref_iter = outer_iter_ref ;
   ss->max_inner_iter = inner_iter_ref ;
   ss->diag_perturb=diag_perturb;*/
   zero_thresh = diag_perturb;
   //status = analyze_kkt();
   //ss->symbolic_analysis();
   k_size = H->ncol + B->nrow + A->nrow;
   skkt_col = k_size;
   sKKTrhs = new double[k_size]();
   kkt_solution = new double[k_size]();
   //diag_perturb = pow(10, -6);

   descent = kkt_solution;
   lagrange_mult_eq = kkt_solution + H->ncol;
   lagrange_mult = kkt_solution + H->ncol + A->nrow;
   //sKKT = ss->A_ord;

   return status;
  }

   int initialize_x_cholmod(){
    for (int i = 0; i < H->ncol; ++i) {
     sKKTrhs[i] = -q[i];
    }
    for (int i = H->ncol, j=0; i < H->ncol+A->nrow; ++i, ++j) {
     sKKTrhs[i] = a[j];
    }
    for (int i = H->ncol+A->nrow; i < k_size; ++i) {
     sKKTrhs[i] = .0;
    }
    //ss->ldl_update_variant =2;
    //print_vec("\n RHS: \n",0,sKKT->ncol,sKKTrhs);
    if(A->nrow == 0)
     L = initial_solve_1(kkt,H,NULL,B->nrow,sKKTrhs,kkt_solution,
                         cm,diag_perturb); //unconstrained solution
    else
     L = initial_solve_1(kkt,H,A,B->nrow,sKKTrhs,kkt_solution,
                         cm,diag_perturb); //unconstrained solution

    //ss->ldl_update_variant =2;
    for (int i = 0; i < H->ncol; ++i) {
     primal_vars[i] = kkt_solution[i];
    }
    for (int i = 0, j = H->ncol; i < A->nrow; ++i, ++j) {
     dual_vars_eq[i] = kkt_solution[j];
    }
  /*  for (int k = 15293; k < 16465; ++k) {
    primal_vars[k] = 0;
   }*/
    //print_vec("\nX0: \n",0,H->ncol,primal_vars);
    //print_vec("\nD1: \n",H->ncol,H->ncol+A->nrow,kkt_solution);
    //print_vec("\nD2: \n",H->ncol+A->nrow,H->ncol+A->nrow+B->nrow,kkt_solution);
    return 1;
   }

   int edit_kkt_cholmod(int add_del, int nxt_col ){
    int row_of_kkt = H->ncol+A->nrow+nxt_col;
    edit_l_factor(row_of_kkt,nxt_col,L,BT,add_del,cm,diag_perturb);
   }

   int solve_kkt_cholmod(solve_type s_type){
    int ret_val = 0;
    chol_solve_updated_kkt(sKKTrhs, L, cm, kkt_solution,
                           kkt_updated,H,B,A,active_set,inner_iter_ref);
    return ret_val;
   }
#endif

/*
 * Creating symbolic info for QP.
 * building super KKT within the solver body.
 */
  int symbolic_QP();

  void reset_symbolic_info();

  void auto_perturbation(size_t total_nnz, double H_norm2, size_t dim,
                         double &reg_diag,
                         int &outer_iter_ref, int &inner_iter_ref);

  void auto_perturbation2(size_t dim, size_t total_nnz, double H_norm2, double range,
                          int n_ineq, double ineq_ratio, double acc_thresh,
                          double &reg_diag,
                          int &outer_iter_ref, int &inner_iter_ref,
                          double &tol);


  void determine_iterations(size_t total_nnz, double acc_thresh,
                            int &outer_iter_ref, int &inner_iter_ref);

  void allocate_workspace();

/*
 * warm-start initilization
 */
  int set_warm_start(double *x0_in, double *dual_in,
                     const std::vector<int> &actve_set,
                     double *dual_in_eq = NULL);

/*
 * warm-start initilization
 */
  int set_warm_start_primal(double *x0_in);

/*
 * Initializes QP with its unconstrained solution.
 */
  int initialize_x();

/*
 * Doing the numerical part, actual iterations.
 */
  int numeric_QP();


/*
 * Update dual vars with descent
 */
  void update_dual(double step);

/*
* Update dual vars with descent
*/
  void update_primal(double step);

/*
 * Update active-set, either dropping or adding a constraint
 */
  int update_active_set(int add_drop, int c_no);

/*
 * The main entry! it solves the QP.
 */
  int solve();

  nasoq_status check_solve_status(int st);

  void detect_solver_name();

/*
* The main entry! it solves the QP using warm-start.
*/
  int solve(double *primal_in, double *dual_in, const std::vector<int> &as);

/*
* The main entry! it solves the QP using primal warm-start.
*/
  int solve(double *primal_in);


/*
 * Returns index number of most violated constraints
 * otherwise, it returns -1 that means the solution is feasible.
 */
  int primal_feasibility(int &nxt_add);

/*
* If return value positive, it is infeasible
* if it is zero, it is feasible.
*/
  int dual_feasibility();

/*
 * Computes dual step length
 */
  double dual_step_length(int &nxt_drop);

/*
 * Computes primal step length
 */
  double primal_step_length(int nxt_active);

  int update_rhs(int nxt_active);

  int solve_kkt_from_scratch();

  int solve_kkt(solve_type s_type);

/*
 * Editing the nxt_col in the KKT matrix.
 */
  void edit_kkt(int add_del, int nxt_col);

/*
 * Computes primal objective
 */
  double compute_primal_obj();

  double quad_form(const CSC *P, const double *x);


/*
 * Computes dual objective
 */
  double compute_dual_obj();

/*
* Computes ||dual_vars * b-Bx||inf
* ||dual_vars_eq * a-Ax||inf
*/
  double complementarity_norm();


/*
 *
 */
  double compute_objective();

/*
 * Computes the norm of Hx + q + B^T*dual_vars
 */
  double lagrangian_residual_norm();

/*
* Computes min(b-Bx,0)
*/
  double constraint_sat_norm();

/*
 * a_i = x_i;
 * b_i = [C*A^{-1}*C^T*\lambda - C*A^{-1}*b + d]_i
 */
  double dual_FB_norm();

/*
 * computes the norm of primal Fischer-Burmeister
 */
  double primal_FB_norm();

/*
 * min(0, \lambda_i)
 */
  double non_negativity_norm();


/*
 * Print QP parametes
 */
  void print(int vb_level = 0);

/*
 * Print QP parameters
 */
  void export_to_file(std::string p_name, std::string vec_head = "");

/*
 * Print log for CSV format
 */
  void print_log();

/*
 *
 */
  void print_qp_range(int verbose = 0);

  int scale_data_2();

  int unscale_data_2();

  int unscale_solution_2();

/*
 *
 */
  int scale_data();

/*
 *
 */
  int unscale_data();

/*
*
*/
  int unscale_solution();


/*
* Doing the numerical part in batch, actual iterations. //FIXME not done
*/
  int numeric_QP_batch();

/*
* Editing the nxt_col in the KKT matrix.
*/
  int edit_kkt_batch(int add_del, std::vector<idx_val> nxt_list);

/*
* Computes primal step length in batch
*/
  double primal_step_length_batch(const std::vector<idx_val> &nxt_list);


/*
* Returns index number of most violated constraints
* otherwise, it returns -1 that means the solution is feasible.
*/
  int primal_feasibility_batch(std::vector<idx_val> &pair_list);

/*
* find the most negative lagrange multiplier
*/
  int dual_feasibility_inner(int &most_neg);

/*
* Update active-set, either dropping or adding a constraint in batch
*/
  int update_active_set_batch(int add_drop,
                              const std::vector<idx_val> nxt_lst);

 };
}
#endif //PROJECT_QP_DUAL_H
