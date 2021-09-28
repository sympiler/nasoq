//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_UTIL_H
#define CHOLOPENMP_UTIL_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>


#include "nasoq/common/def.h"
#include "nasoq/common/transpose_unsym.h"

namespace nasoq {
/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
 bool readMatrix_old(std::string fName, int &n, int &NNZ, int *&col,
                     int *&row, double *&val);

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
 bool readMatrix(std::string fName, size_t &n, size_t &NNZ, int *&col,//FIXME change col type to size_t
                 int *&row, double *&val);


/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
 bool readMatrix_rect(std::string fName, size_t &n_row, size_t &n_col,
                      size_t &NNZ, int *&col,//FIXME change col type to size_t
                      int *&row, double *&val);

/*
 *
 */
 void read_dense(std::string fname, int &n_row, int &n_col, double *&mat);

/*
 *
 */
 int write_dense(std::string fname, int n_row, int n_col, double *mat);

/*
 * expand the empty columns by one zero
 * e.g., in indefinite matrices where empty columns exist
 * returns false if the matrix is not expanded.
 */
 bool expandMatrix(size_t n, size_t NNZ, const int *col,
                   const int *row, const double *val,
                   size_t &newNNZ, int *&newCol,//FIXME change col type to size_t
                   int *&newRow, double *&newVal,
                   double insDiag = 0);

/*
 * reading a ordering from PasTiX stored format.
 */
 int read_vector(std::string fName, size_t n, double *perm);

/*
 * writing a vector to a file.
 */
 int write_vector(std::string fName, size_t n,
                  double *vec_vals, std::string header = "",
                  int prec = 30);

/*
 * reading a ordering from PasTiX stored format.
 */
 bool readOrdering(std::string fName, size_t n, size_t *perm);

 bool enableColdCache(int n, std::ifstream &f);

/*
 * For linear solver using A and A'
 */

 void rhsInit_linearSolver(int n, int *Ap, int *Ai, double *Ax, //A
                           int *Bp, int *Bi, double *Bx, //AT
                           double *b);

/*
 * Assuming we have only lower part of a symmetric matrix
 * return false if rhs is zero
 */
 bool rhsInit_linearSolver_onlylower(int n, int *Ap, int *Ai, double *Ax, //A
                                     double *b);

/*
 * For triangular solve
 */

 void rhsInit(int n, int *Ap, int *Ai, double *Ax, double *b);

/*
 * RHS initilization for blocked
 */

 void rhsInitBlocked(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, double *Ax, double *b);

/*
 * RHS initilization for blocked L'
 */

 void rhsInitBlockedLT(size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP, double *Ax, double *b);

/*
 * Making full symmetric matrix from lower and upper
 */

 void make_full(int ncol, int nnz, int *Ap, int *Ai, double *Ax, //A
                int *ATransp, int *ATransi, double *ATransx,
                int &nnzFull, int *&AFullp, int *&AFulli, double *&AFullx);

/*
 * Testing lower triangular solve to make sure all x elements are ONE.
 * Warning: This works only for when RHS does not have zero
 */

 int testTriangular(size_t n, const double *x, double epsilon = 1e-3);



/*
 * converting CSC to CSR
 */  //TODO
/*
 void csc2csr(int n, int *Ap, int *Ar, double *Ax,
              int *row_ptr, int *col_idx, double *Bx) {
  int *row_len = new int[n]();
  int *finger = new int[n]();
  for (int i = 0; i < Ap[i + 1]; ++i) {
   row_len[Ar[i]]++;
  }
  col_idx[0] = 0;
  for (int l = 1; l < n + 1; ++l) {
   col_idx[l] = col_idx[l - 1] + row_len[l];
  }
  for (int k = 0; k < n; ++k) {
   for (int i = Ap[i]; i < Ap[i + 1]; ++i) {
    //Bx[col_idx[]]
   }
  }
 }
*/


/*
 * Converting BCSC to CSC
 */
 int bcsc2csc(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col, double *Ax,
   //out
   int *Cp, int *Ci, double *Cx
 );


/*
 * Converting BCSC to CSC, removing all zeros, assuming
 * fill-ins are initially zero.
 * Warning: could be problematic in some exceptional cases where,
 * exact numerical cancellation happens.
 */
 int bcsc2csc_aggressive(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col, double *Ax,
   //out
   //MKL_INT *Cp, MKL_INT *Ci, double *Cx
   size_t *Cp, int *Ci, double *Cx
 );


 int check_row_idx_l(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col
 );

 int bcsc2csc_aggressive_int(
   //in
   size_t n, size_t nBlocks, size_t *Ap, int *Ai, size_t *AiP,
   int *sup2col, double *Ax,
   //out
   //MKL_INT *Cp, MKL_INT *Ci, double *Cx
   int *Cp, int *Ci, double *Cx
 );

 void swapWSet(std::vector<std::vector<int>> &List, int i, int j);


/*
 * Check whether the matrix is stored in full
 * return the lower half if it is full
 * Makes sense only for symmetric matrices.
 */
 CSC *computeLowerTriangular(CSC *A);


 void print_csc(std::string beg, size_t n, int *Ap, int *Ai, double *Ax);

 void print_csc_l(std::string beg, size_t n,
                  long *Ap, int *Ai, double *Ax);

 template<class type>
 void print_vec(std::string header, int beg,
                int end, type *vec) {

  std::cout << header;
  for (int i = beg; i < end; ++i) {
   std::cout << std::setprecision(15) << vec[i] << ", ";
  }
  std::cout << "\n";
 }

/*
 * Print dense matrix, stored col-wise
 *
 */
 template<class type>
 void print_mat(std::string header, int row_no,
                int col_no, type *mat) {
  std::cout << header;
  for (int i = 0; i < col_no; ++i) {
   for (int j = 0; j < row_no; ++j) {
    std::cout << std::setprecision(15) << mat[i * row_no + j] << ", ";
   }
   std::cout << "\n";
  }
  std::cout << "\n";
 }

 int *compute_inv_perm(int n, int *perm, int *pinv);


 void combine_perms(int n, int *perm1, int *perm2, int *result);

 int zero_diag_removal(size_t n, size_t NNZ, const int *col,
                       const int *row, double *val,
                       double threshold, double reg_diag);

/*
 * return 0 if equal
 */
 int is_vector_equal(int n, double *v1, double *v2,
                     double eps = std::numeric_limits<double>::min());

/*
 * if all vector values are zero
 */
 int is_vector_zero(int n, double *v, double threshold);


/*
 * returns True:1 if equal
 */
 int is_value_equal(double a, double b, double eps);


/*
 * removes duplicates of a vector
 */
 void make_unique(std::vector<int> &sn_arrays);

/*
 * Gathering rows from a CSC matrix
 */
 void gather_rows(std::vector<int> row_ids,
                  size_t An, size_t Am,
                  int *Ap,
                  int *Ai, double *Ax,
                  size_t &Bn, size_t &Bm, size_t &Bnz, int *&Bp,
                  int *&Bi, double *&Bx);


/*
* Build super matrix from A and B,  All are stored in CSC.
*/
 int build_super_matrix_test(CSC *A, CSC *B, double reg_diag,
                             size_t &SM_size, size_t &SM_nz,
                             int *&SMp, int *&SMi, double *&SMx,
                             double *&sm_rhs);

 int make_lower(size_t An, size_t Anz, int *Ap, int *Ai, double *Ax,
                size_t &Bn, size_t &Bnz, int *&Bp, int *&Bi, double *&Bx);

 void compressed_set_to_vector(int set_size, int *set_p, int *set_n,
                               std::vector<std::vector<int>> &res_set);

 void max_min_spmat(size_t An, int *Ap, int *Ai, double *Ax,
                    double &max_val, double &min_val);


 bool parse_args(int argc, const char *argv[], std::map<std::string, std::string> &qp_args);

/*
 bool parse_args(int argc, char **argv, std::map<std::string, std::string> &qp_args) {
  const char *const short_opts = "m:o:l:a:b:c:d:n:p:r:e:t:h";
  const option long_opts[] = {
    {"mode",       required_argument, nullptr, 'm'},
    {"objective",  required_argument, nullptr, 'o'},
    {"linear",     required_argument, nullptr, 'l'},
    {"eq",         required_argument, nullptr, 'a'},
    {"eqb",        required_argument, nullptr, 'b'},
    {"ineq",       required_argument, nullptr, 'c'},
    {"ineqb",      required_argument, nullptr, 'd'},
    {"nasoq",      required_argument, nullptr, 'n'},
    {"perturb",    required_argument, nullptr, 'p'},
    {"refinement", required_argument, nullptr, 'r'},
    {"epsilon",    required_argument, nullptr, 'e'},
    {"toli",       required_argument, nullptr, 't'},
    {"help",       no_argument,       nullptr, 'h'},
    {nullptr,      no_argument,       nullptr, 0}
  };

  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  while (true) {
   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
   if (-1 == opt)
    break;

   switch (opt) {
    case 'm':
     qp_args.insert(std::pair<std::string, std::string>("mode", optarg));
     break;
    case 'o':
     qp_args.insert(std::pair<std::string, std::string>("H", optarg));
     break;
    case 'l':
     qp_args.insert(std::pair<std::string, std::string>("q", optarg));
     break;
    case 'a':
     qp_args.insert(std::pair<std::string, std::string>("A", optarg));
     break;
    case 'b':
     qp_args.insert(std::pair<std::string, std::string>("b", optarg));
     break;
    case 'c':
     qp_args.insert(std::pair<std::string, std::string>("C", optarg));
     break;
    case 'd':
     qp_args.insert(std::pair<std::string, std::string>("d", optarg));
     break;
    case 'n':
     qp_args.insert(std::pair<std::string, std::string>("nasoq", optarg));
     break;
    case 'p':
     qp_args.insert(std::pair<std::string, std::string>("perturbation", optarg));
     break;
    case 'r':
     qp_args.insert(std::pair<std::string, std::string>("iterations", optarg));
     break;
    case 'e':
     qp_args.insert(std::pair<std::string, std::string>("epsilon", optarg));
     break;
    case 't':
     qp_args.insert(std::pair<std::string, std::string>("tolerance", optarg));
     break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
     print_help();
     break;
   }
  }
  return true;
 }
*/

/*
 bool parse_nasoq_args(int argc, char **argv, std::map<std::string, std::string> &qp_args) {
  const char *const short_opts = "i:o:l:d:v:p:r:e:t:h";
  const option long_opts[] = {
    {"input",      required_argument, nullptr, 'i'},
    {"output",      required_argument, nullptr, 'o'},
    {"log",      required_argument, nullptr, 'l'},
    {"header",      required_argument, nullptr, 'd'},
    {"variant",      required_argument, nullptr, 'v'},
    {"perturb",    required_argument, nullptr, 'p'},
    {"refinement", required_argument, nullptr, 'r'},
    {"epsilon",    required_argument, nullptr, 'e'},
    {"toli",       required_argument, nullptr, 't'},
    {"help",       no_argument,       nullptr, 'h'},
    {nullptr,      no_argument,       nullptr, 0}
  };
  auto print_help = [](){std::cout<<"Input argument is wrong, you need at least an input QP file ";};
  while (true) {
   const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

   if (-1 == opt)
    break;
   switch (opt) {
    case 'i':
     qp_args.insert(std::pair<std::string, std::string>("input", optarg));
     break;
    case 'o':
     qp_args.insert(std::pair<std::string, std::string>("output", optarg));
     break;
     case 'l':
     qp_args.insert(std::pair<std::string, std::string>("log", optarg));
     break;
    case 'd':
     qp_args.insert(std::pair<std::string, std::string>("header", optarg));
     break;
     case 'v':
     qp_args.insert(std::pair<std::string, std::string>("variant", optarg));
     break;
    case 'p':
     qp_args.insert(std::pair<std::string, std::string>("perturbation", optarg));
     break;
    case 'r':
     qp_args.insert(std::pair<std::string, std::string>("iterations", optarg));
     break;
    case 'e':
     qp_args.insert(std::pair<std::string, std::string>("epsilon", optarg));
     break;
    case 't':
     qp_args.insert(std::pair<std::string, std::string>("tolerance", optarg));
     break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
     print_help();
     break;
   }
  }
  if(qp_args.find("input") == qp_args.end()){
   print_help();
   return false;
  }
  return true;
 }
*/



}
#endif //CHOLOPENMP_UTIL_H
