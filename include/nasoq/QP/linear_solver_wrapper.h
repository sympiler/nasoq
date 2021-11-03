//
// Created by kazem on 3/5/19.
//

#ifndef PROJECT_LINEAR_SOLVER_WRAPPER_H
#define PROJECT_LINEAR_SOLVER_WRAPPER_H

#include <chrono>
#include <vector>

#include "nasoq/common/def.h"

//#undef SYM_REMOV

namespace nasoq {
 enum solve_type {
  FULL = 0, REFACTOR = 1, UPDATE = 2, SOLVE = 3
 };

 struct profiling_solver_info {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;

  double fact_time, analysis_time, solve_time, iter_time, ordering_time;
  double update_time;
  double piv_reord;
  double *timing_chol;

  profiling_solver_info(int nt);

  profiling_solver_info(int nt, double ft, double at, double st,
                        double it, double ot, double ut, double pr);

  ~profiling_solver_info();

  std::chrono::time_point<std::chrono::system_clock> tic();

  std::chrono::time_point<std::chrono::system_clock> toc();

  double elapsed_time(std::chrono::time_point<std::chrono::system_clock> beg,
                      std::chrono::time_point<std::chrono::system_clock> lst);

  void print_profiling();


 };

 struct perturbed_value {
  double per_value;
  int col_idx;
 };

 struct SolverSettings {
  int solver_mode; //0 is normal solve, 1 is row/col addition
  BCSC *L;
  double *valL, *d_val, *x, *rhs;
  CSC *A, *SM; // SuperMatrix
  CSC *AorSM; // Pointer to A or SM
  CSC *A_ord, *AT_ord; //reordered matrix, either SM or A
  CSC *B, *BT;
  CSC *C, *CT; // These are for QP mode.
  size_t dim1;
  size_t base; //base is zero for mode 0 and SM->ncol - B->nrow for mode 1
  int *perm_piv, *pinv;
  int a_consistent; // Is A ordered according to L?
  int ldl_variant, ldl_update_variant;
  double *sm_rhs, *sm_solution;
  int *l_pb, *l_pe, *l_i;//simplicial format
  double *l_x;
  int is_super, simplicial_alloc;
  size_t ws_int_size;
  size_t ws_dbl_size;

  //parallel variables
  int n_level, n_par;
  int *level_ptr, *level_set, *par_ptr, *par_set;
  int n_level_s, n_par_s;
  int *level_ptr_s, *par_ptr_s, *par_set_s;
  int *s_level_ptr, *s_level_set;
  int s_level_no;
  int *prune_ptr, *prune_set;
  int *extra_cols;
  //sparsity related
  int n_relax[3];
  double z_relax[3];
  int max_sup_wid, max_col, status=0;
  int ordering_type; //TODO to be used later.
  //etree info
  int *atree, *etree, *etree_mod;
  int *visible_cnt; //number of invisible col of a supernode
  int *child_ptr, *child_no, *num_child;
  int *child_sn_ptr, *child_sn_no, *num_sn_child;
  std::vector<std::vector<int>> children_vec;
  std::vector<int> detached_nodes;
  int col_del, to_del;
  //Numeric related
  double reg_diag;
  int num_pivot;
  // For iterative refinement
  int max_iter, max_inner_iter;
  int num_ref_iter; //output
  int req_ref_iter;
  double tol_abs, tol_rel;
  std::vector<perturbed_value> perturbed_diags;
  int regularization;

  //architecture info
  int cost_param, level_param, final_seq_node;
  int chunk, num_thread, thread_thresh;

  //Profiling info FIXME: define profiling class
  profiling_solver_info *psi;

  //Update/Downdate
  std::vector<int> modified_sns;
  bool *marked, *visible_sn;
  double *ws; // workspace for double computations
  double *ws_zeroed; //after each process, this is zero
  int *ws_int; // workspace for integer computations
  double *extra_rhs;
  //Internal vars
  int remove_trans;

  //Norms
  double A_l1, x_l1, rhs_l1, res_l1;
  double bwd_err;

  SolverSettings(CSC *Amat, double *rhs_in);

  SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat);

  SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, double *b);

  SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, CSC *BTmat);

  SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, CSC *BTmat,
                 CSC *Cmat, CSC *CTmat);

  ~SolverSettings();

  /*
   * Default setting of solver
   */
  void default_setting();

  /*
  * Build superKKT from transpose of Hessian and
  * the constraints. All stored in CSC.
  */
  int build_super_matrix();

  /*
   *
   */
  void find_perturbation(double tol);

  void apply_perturbation(double tol);

  /*
   * Applies perturbation to the reorded system
   */
  void add_perturbation(double tol);

  void remove_perturbation(double tol);


  /*
   *
   */
  int symbolic_analysis();

  /*
   *Assuming objective and constraint matrices are fixed.
   */
  void reset_symbolic_factor();

  /*
   *
   */
  void architecture_related_params();

/*
 * Allocated all intermediate required memory
 */
  void allocate_workspace();

  /*
   *
   */
  int numerical_factorization();

  /*
  * deleting node k from atree assuming the col is in supernode
   * with size of 1.
  */
  int delete_node_tree_simple(int k);

  /*
 * Adding node k to atree
  * We need to access original supermatrix etree
   * This routine assumes the added supernode has only one col
 */
  int add_node_tree_simple(int k);

  void bfs_tree(int cur_sn);

  /*
   * deleting node k from atree
   */
  int delete_node_tree(int k);

  /*
  * Adding node k to atree
   * We need to access original supermatrix etree
  */
  int add_node_tree(int k);


/*
 * This is an expert routine, used for QP
 * Edit the list of col/rows of matrix B
 * in add_drop_const.
 * rhs_values will change the corresponding values of the RHS
 */
  int add_del_matrix_qp(int add_del, std::vector<int> add_drp_const);

  /*
   * Updates the linear system by add rows/cols in add_drp_const vector
   * And also udpate corresponding location by the given input rhs
   */
  int update_somod(std::vector<int> add_drp_const, const int n_rhs, const
  double *new_rhs);

  /*
   * Removes (makes zero) rows/cols in add_drp_const
   */
  int downdate_somod(std::vector<int> add_drp_const, int n_rhs);

  /*
   * Editing an existing row
   */
  int edit_matrix();

  /*
   *
   */
  int edit_rhs();

  /*
   *
   */
  double *edit_solve_rhs();

  /*
  * only updates L, assuming matrix A is already edited and
   * marked or modified_sn are already set.
  */
  int update_factorization();

/*
 * Edits the row/col in original symmetric matrix A and updates L
 */
  int edit_update_factorization(int add_del, std::vector<int> add_drp_const);

  double *update_solve();


  /*
   * Use the computed factors and solve the system for an internal RHS for a given rhs_in
   */
  double *solve_only();
  double *solve_only(const int n_rhs, double *rhs_in);

  double *iterative_ref_only();

  void reorder_matrix();

  /*
   * Computes all required norms
   */
  void compute_norms();

  double backward_error();

/*
  * creating non-supernodal factor by preserving the
  * value array intact.
  */
  void convert_supernode_to_simplicial();

/*
 * Copy values from L supernodal to L simplicial
 */
  void copy_lsuper_to_l();

  /*
   * Set diagonals, used for perturbation on reoredered matrix
   * set it to zero if d == 0 otherwise use the value in the
   * perturbed_diags array
   */
  void set_diags(double d);


  int check_ldlt_factor();

 };

}
#endif //PROJECT_LINEAR_SOLVER_WRAPPER_H
