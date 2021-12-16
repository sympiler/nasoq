//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/QP/updown_test.h"

#ifdef MKL_BLAS
#include "mkl.h"
#endif
#include "nasoq/common/Util.h"
#include "nasoq/common/Transpose.h"
#include "nasoq/QP/linear_solver_wrapper.h"

namespace nasoq {

 int build_super_matrix_with_eq(CSC *A, CSC *B, CSC *C, double reg_diag, size_t &SM_size, size_t &SM_nz, int *&SMp,
                                int *&SMi, double *&SMx, double *&sm_rhs) {
  int status=0;
  if (B->nrow > 0) {
   if (C == NULL || C->nrow == 0) {
    SM_size = A->ncol + B->nrow;
    SMp = new int[SM_size + 1];
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (B->p[i] - B->p[i - 1]);
    }
   } else {
    SM_size = A->ncol + B->nrow + C->nrow;
    SMp = new int[SM_size + 1];
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (B->p[i] - B->p[i - 1]) +
              (C->p[i] - C->p[i - 1]);
    }
   }
  } else {
   if (C == NULL || C->nrow == 0) {
    SM_size = A->ncol;
    SMp = new int[SM_size + 1];
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]);
    }
   } else {
    SM_size = A->ncol + C->nrow;
    SMp = new int[SM_size + 1];
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (C->p[i] - C->p[i - 1]);
    }
   }
  }
  //Adding diagonal for columns with zero values.
  for (int k = A->ncol + 1; k < SM_size + 1; ++k) {
   SMp[k] = SMp[k - 1] + 1;
  }
  SM_nz = SMp[SM_size];
  SMi = new int[SM_nz];
  SMx = new double[SM_nz]();

  int base1 = A->ncol;
  size_t stp = 0;
  for (int j = 0; j < A->ncol; ++j) {
   stp = SMp[j];
   //Adding Hessian
   for (int i = A->p[j]; i < A->p[j + 1]; ++i) {
    SMi[stp] = A->i[i];
    SMx[stp] = A->x[i];
    stp++;
   }
   base1 = A->ncol;
   //Adding equalities
   if (C != NULL || C->nrow > 0) {
    if (C->nrow > 0) {
     for (int i = C->p[j]; i < C->p[j + 1]; ++i) {
      SMi[stp] = base1 + C->i[i];
      SMx[stp] = C->x[i];
      stp++;
     }
    }
    base1 = A->ncol + C->nrow;
   }
   //Adding inequalities
   if (B->nrow > 0) {
    for (int i = B->p[j]; i < B->p[j + 1]; ++i) {
     SMi[stp] = base1 + B->i[i];
     SMx[stp] = B->x[i];
     stp++;
    }
   }
   assert(stp == SMp[j + 1]);
  }
  //Putting a small value in diagonals
  base1 = A->ncol;
  for (int l = SMp[A->ncol], j = 0; l < SM_nz; ++l, ++j) {
   //SMx[l] = diag_perturb;
   SMi[l] = base1 + j;
  }

  sm_rhs = new double[SM_size]();
  //sm_solution = new double[SM->ncol]();
//  CSC *sKKTt = ptranspose(SM,2,NULL,NULL,0,status);
  //print_csc("skkt\n",SM_size,SMp,SMi,SMx);
  //print_csc("\nskkt Trans\n",sKKTt->ncol,sKKTt->p,sKKTt->i,sKKTt->x);
  //std::cout<<"\n";
  return status;
 }

 void build_super_solve_with_eq(CSC *H, const CSC *A, CSC *B, const double *rhs, double reg_diag,
                                const std::vector<int> &mod_col, double *x_all, int outer_iter, int inner_iter,
                                double stop_tol) {
  CSC *SM_test;
  double *rhs_test;
  SolverSettings *stmp_all;
  CSC *A_red = new CSC;
  SM_test = new CSC;
  gather_rows(mod_col, A->nrow, A->ncol, A->p, A->i, A->x,
              A_red->nrow, A_red->ncol, A_red->nzmax,
              A_red->p, A_red->i, A_red->x);
  build_super_matrix_with_eq(H, A_red, B, 0, SM_test->ncol,
                             SM_test->nzmax,
                             SM_test->p, SM_test->i, SM_test->x,
                             rhs_test);

  for (int l = 0; l < H->ncol + A_red->nrow + B->nrow; ++l) {
   rhs_test[l] = rhs[l];
  }

  SM_test->nrow = SM_test->ncol;
  SM_test->stype = -1;
  SM_test->xtype = CHOLMOD_REAL;
  SM_test->packed = TRUE;
  SM_test->sorted = TRUE;
  stmp_all = new SolverSettings(SM_test, rhs_test);

  //print_csc("full :\n",SM_test->ncol,SM_test->p,SM_test->i,SM_test->x);
  SM_test->nrow = SM_test->ncol;
  SM_test->stype = -1;
  SM_test->xtype = CHOLMOD_REAL;
  SM_test->packed = TRUE;
  SM_test->sorted = TRUE;

  stmp_all->level_param = 0;
  stmp_all->final_seq_node = 4;
  stmp_all->ldl_variant = 1;
  stmp_all->req_ref_iter = outer_iter;
  stmp_all->max_inner_iter = inner_iter;
  stmp_all->reg_diag = reg_diag;
  stmp_all->tol_rel = stmp_all->tol_abs = stop_tol;
  stmp_all->solver_mode = 0;
  stmp_all->symbolic_analysis();
  stmp_all->numerical_factorization();
  stmp_all->solve_only();
  //stmp_all->backward_error();
  for (int k = 0; k < SM_test->ncol; ++k) {
   x_all[k] = stmp_all->x[k];
  }
  //std::cout<< " BWD: "<<stmp_all->backward_error()<<"\n";
  //stmp_all->psi->print_profiling();


  delete[]stmp_all->x;
  delete stmp_all;
  allocateAC(A_red, 0, 0, 0, FALSE);
  allocateAC(SM_test, 0, 0, 0, FALSE);
  delete[]rhs_test;
 }

#ifdef MKL_BLAS
 void build_super_solve_with_eq_mkl(CSC *H, const CSC *A, CSC *B, const double *rhs, double reg_diag,
                                    const std::vector<int> &mod_col, double *x_all, int outer_iter, int inner_iter,
                                    double stop_tol) {

  CSC *SM_test;
  double *rhs_test;
  SolverSettings *stmp_all;
  CSC *A_red = new CSC;
  SM_test = new CSC;
  double *a;
  int *ja;
  int *ia;
  gather_rows(mod_col, A->nrow, A->ncol, A->p, A->i, A->x,
              A_red->nrow, A_red->ncol, A_red->nzmax,
              A_red->p, A_red->i, A_red->x);
  build_super_matrix_with_eq(H, A_red, B, 0, SM_test->ncol,
                             SM_test->nzmax,
                             SM_test->p, SM_test->i, SM_test->x,
                             rhs_test);

  for (int l = 0; l < H->ncol + A_red->nrow + B->nrow; ++l) {
   rhs_test[l] = rhs[l];
  }

  SM_test->nrow = SM_test->ncol;
  SM_test->stype = -1;
  SM_test->xtype = CHOLMOD_REAL;
  SM_test->packed = TRUE;
  SM_test->sorted = TRUE;
  //print_csc("full :\n",SM_test->ncol,SM_test->p,SM_test->i,SM_test->x);
  SM_test->nrow = SM_test->ncol;
  SM_test->stype = -1;
  SM_test->xtype = CHOLMOD_REAL;
  SM_test->packed = TRUE;
  SM_test->sorted = TRUE;
  SM_test->nzmax = SM_test->p[SM_test->ncol];
  int status = 0;
  CSC *SMT = ptranspose(SM_test, 2, NULL, NULL, 0, status);
  CSC *H_full = new CSC;
  int nnz_full = 0;
  make_full(SM_test->ncol, SM_test->nzmax, SM_test->p, SM_test->i,
            SM_test->x, SMT->p, SMT->i, SMT->x, nnz_full,
            H_full->p, H_full->i, H_full->x);
  H_full->nrow = H_full->ncol = SM_test->ncol;
  H_full->stype = 0;
  H_full->packed = 1;
  H_full->nzmax = H_full->p[H_full->ncol];

  a = H_full->x;
  ia = H_full->p;
  ja = H_full->i;
  int n = SM_test->ncol;

  MKL_INT mtype = 11;       /* Real unsymmetric matrix */
  /* RHS and solution vectors. */
  //double b[8], x[8];
  double *b = rhs_test;
  double *x = (double *) mkl_calloc(n, sizeof(double), 64);
  MKL_INT nrhs = 1;     /* Number of right hand sides. */
  long int pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  MKL_INT i;
  double ddum;          /* Double dummy */
  MKL_INT idum;         /* Integer dummy. */


/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
   iparm[i] = 0;
  }
  iparm[0] = 1;         /* No solver default */
  iparm[1] = 2;         /* Fill-in reordering from METIS */
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
// iparm[4] = 1;         /* user fill-in reducing permutation is used*/
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 4;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 9;        /* Perturb the pivot elements with 1E-8 */
  iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* A^TX=B */
  iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = 1;       /* Output: Mflops for LU factorization */
  iparm[19] = inner_iter;        /* Output: Numbers of CG Iterations */
  iparm[20] = 1; /*using bunch kaufman pivoting*/
  iparm[55] = 0; /*Diagonal and pivoting control., default is zero*/

  //
  iparm[26] = 1;
  //iparm[23] = 1; //TODO: Later enable to se if the parallelism is better
  iparm[34] = 1;
  //Because iparm[4]==0 so:
  iparm[30] = 0;
  iparm[35] = 0;
  maxfct = 1;           /* Maximum number of numerical factorizations. */
  mnum = 1;         /* Which factorization to use. */
  msglvl = 0;           /* Print statistical information in file */
  error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
   pt[i] = 0;
  }

/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */

  phase = 11;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
   printf("\nERROR during symbolic factorization: %d", error);
   exit(1);
  }

/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
   printf("\nERROR during numerical factorization: %d", error);
   exit(2);
  }
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
  phase = 33;
  iparm[7] = 4;         /* Max numbers of iterative refinement steps. */

  PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
          &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, b, x, &error);
  if (error != 0) {
   printf("\nERROR during solution: %d", error);
   exit(3);
  }


  for (int k = 0; k < SM_test->ncol; ++k) {
   x_all[k] = x[k];
  }
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
  phase = -1;           /* Release internal memory. */
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
          &n, &ddum, ia, ja, NULL, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error);

  allocateAC(A_red, 0, 0, 0, FALSE);
  allocateAC(SM_test, 0, 0, 0, FALSE);
  allocateAC(H_full, 0, 0, 0, FALSE);
  allocateAC(SMT, 0, 0, 0, FALSE);
  delete[]rhs_test;
  mkl_free(x);
 }

#endif
}