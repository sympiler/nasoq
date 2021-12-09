//
// Created by kazem on 2019-12-05.
//

#ifndef PARS_CHOLMOD_UTILS_H
#define PARS_CHOLMOD_UTILS_H
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cholmod.h>
#include <cholmod_function.h>
#include <omp.h>
#include "nasoq/common/def.h"
#include "nasoq/QP/updown_test.h"

namespace nasoq {
 void print_sparse_l(cholmod_factor *L) {
  long *Lp = (long *) L->p;
  long *Li = (long *) L->i;
  double *Lx = (double *) L->x;
  for (int i = 0; i < L->n; ++i) {
   for (int j = Lp[i]; j < Lp[i + 1]; ++j) {
    std::cout << i << " " << Li[j] << " " << Lx[j] << "\n";
   }
  }
  std::cout << "\n";
 }

/*
* Build super matrix from A and B,  All are stored in CSC.
 * C is eq
 * B is ineq
*/
 int build_super_matrix_with_eq_chol(size_t total_size,
                                     CSC *A, CSC *B, CSC *C, double reg_diag,
                                     size_t &SM_size, size_t &SM_nz,
                                     int *&SMp, int *&SMi, double *&SMx,
                                     double *&sm_rhs) {
  int status;
  SM_size = total_size;
  if (B->nrow > 0) {
   if (C == NULL || C->nrow == 0) {
    //SM_size = A->ncol + B->nrow ;
    SMp = new int[SM_size + 1];
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (B->p[i] - B->p[i - 1]);
    }
   } else {
    //SM_size = A->ncol + B->nrow + C->nrow;
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
    //SM_size = A->ncol ;
    SMp = new int[SM_size + 1];
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]);
    }
   } else {
    //SM_size = A->ncol + C->nrow;
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
   SMx[l] = reg_diag;
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

/******************************************************************************/

 int pmgmres_ldlt_auto_chol(int n, int nz_num, int ia[], int ja[], double a[],
                            cholmod_factor *L, cholmod_common *cm,
                            double x[], double rhs[], int itr_max, int mr,
                            double tol_abs, double tol_rel, int sorted = 1, int blocked = 0)

/******************************************************************************/

 {
  double av;
  double *c;
  double delta = 1.0e-03;
  double *g;
  double **h;
  double htmp;
  int i;
  int itr;
  int itr_used;
  int j;
  int k;
  int k_copy;
  double mu;
  double *r;
  double rho;
  double rho_tol;
  double *s;
  double **v;
  int verbose = 0;
  double *y;
  double alp[2] = {1.0, 0};
  double bet[2] = {0.0, 0};
  itr_used = 0;
  double norm_inf = 0;
  // double workspace: 4*(mr+1) + n + mr(mr-1) + mr(n-1)
  c = (double *) malloc((mr + 1) * sizeof(double));
  g = (double *) malloc((mr + 1) * sizeof(double));
  h = dmatrix(0, mr, 0, mr - 1);
  r = (double *) malloc(n * sizeof(double));
  s = (double *) malloc((mr + 1) * sizeof(double));
  v = dmatrix(0, mr, 0, n - 1);
  y = (double *) malloc((mr + 1) * sizeof(double));
  int xtype = CHOLMOD_REAL;
  cholmod_dense *B_o = cholmod_l_zeros(n, 1, xtype, cm);
  cholmod_dense *Y;
  double *Yx;
  double *B_ox = (double *) B_o->x;

  for (itr = 0; itr < itr_max; itr++) {
   if (!sorted)
    ax_cr(n, nz_num, ia, ja, a, x, r);
   else
    spmv_csc_sym_one_int(n, ia, ja, a, -1, alp, bet, 1, x, r);

   for (i = 0; i < n; i++) {
    r[i] = rhs[i] - r[i];
    B_ox[i] = r[i];
   }
   //norm_inf = norm_dense();
   //print_sparse_l(L);
   //print_vec("dd\n",0,n,B_ox);
   cholmod_dense *TMP = cholmod_l_ones(n, 1, CHOLMOD_REAL, cm);
   Y = cholmod_l_solve(1, L, TMP, cm);
   Yx = (double *) TMP->x;
   for (int m = 0; m < n; ++m) {
    r[m] = Yx[m];
   }
   cholmod_l_free_dense(&Y, cm);
   //print_vec("dd\n",0,n,r);
   rho = sqrt(r8vec_dot(n, r, r));
   if (rho < tol_abs)
    return itr_used;
   if (itr == 0) {
    rho_tol = rho * tol_rel;
   }

   for (i = 0; i < n; i++) {
    v[0][i] = r[i] / rho;
   }
   g[0] = rho;
   for (i = 1; i < mr + 1; i++) {
    g[i] = 0.0;
   }
   for (i = 0; i < mr + 1; i++) {
    for (j = 0; j < mr; j++) {
     h[i][j] = 0.0;
    }
   }
   for (k = 0; k < mr; k++) {
    k_copy = k;
    if (!sorted)
     ax_cr(n, nz_num, ia, ja, a, v[k], v[k + 1]);
    else
     spmv_csc_sym_one_int(n, ia, ja, a, -1, alp, bet, 1, v[k], v[k + 1]);

    for (i = 0; i < n; i++) {
     B_ox[i] = v[k + 1][i];
    }
    Y = cholmod_l_solve(1, L, B_o, cm);
    Yx = (double *) Y->x;
    for (int m = 0; m < n; ++m) {
     v[k + 1][m] = Yx[i];
    }
    cholmod_l_free_dense(&Y, cm);

    av = sqrt(r8vec_dot(n, v[k + 1], v[k + 1]));
    for (j = 0; j <= k; j++) {
     h[j][k] = r8vec_dot(n, v[k + 1], v[j]);
     for (i = 0; i < n; i++) {
      v[k + 1][i] = v[k + 1][i] - h[j][k] * v[j][i];
     }
    }
    h[k + 1][k] = sqrt(r8vec_dot(n, v[k + 1], v[k + 1]));

    if ((av + delta * h[k + 1][k]) == av) {
     for (j = 0; j < k + 1; j++) {
      htmp = r8vec_dot(n, v[k + 1], v[j]);
      h[j][k] = h[j][k] + htmp;
      for (i = 0; i < n; i++) {
       v[k + 1][i] = v[k + 1][i] - htmp * v[j][i];
      }
     }
     h[k + 1][k] = sqrt(r8vec_dot(n, v[k + 1], v[k + 1]));
    }

    if (h[k + 1][k] != 0.0) {
     for (i = 0; i < n; i++) {
      v[k + 1][i] = v[k + 1][i] / h[k + 1][k];
     }
    }

    if (0 < k) {
     for (i = 0; i < k + 2; i++) {
      y[i] = h[i][k];
     }
     for (j = 0; j < k; j++) {
      mult_givens(c[j], s[j], j, y);
     }
     for (i = 0; i < k + 2; i++) {
      h[i][k] = y[i];
     }
    }
    mu = sqrt(h[k][k] * h[k][k] + h[k + 1][k] * h[k + 1][k]);
    if (mu == 0)
     return itr_used;
    c[k] = h[k][k] / mu;
    s[k] = -h[k + 1][k] / mu;
    h[k][k] = c[k] * h[k][k] - s[k] * h[k + 1][k];
    h[k + 1][k] = 0.0;
    mult_givens(c[k], s[k], k, g);

    rho = fabs(g[k + 1]);

    itr_used = itr_used + 1;

    if (rho <= rho_tol && rho <= tol_abs) {
     break;
    }
   }

   k = k_copy;
   if (h[k][k] == 0)
    return itr_used;
   y[k] = g[k] / h[k][k];
   for (i = k - 1; 0 <= i; i--) {
    y[i] = g[i];
    for (j = i + 1; j < k + 1; j++) {
     y[i] = y[i] - h[i][j] * y[j];
    }
    y[i] = y[i] / h[i][i];
   }
   for (i = 0; i < n; i++) {
    for (j = 0; j < k + 1; j++) {
     x[i] = x[i] + v[j][i] * y[j];
    }
   }
   if (rho <= rho_tol && rho <= tol_abs) {
    break;
   }
  }


/*
  Free memory.
*/
  free(c);
  free(g);
  free_dmatrix(h, 0, mr, 0, mr - 1);
  free(r);
  free(s);
  free_dmatrix(v, 0, mr, 0, n - 1);
  free(y);

  return itr_used;
 }


/*
int pmgmres_ldlt_cr_cholmod ( int n, int nz_num, int ia[], int ja[], double a[],
                      cholmod_factor *L, cholmod_common *cm,
                      double x[], double rhs[], int itr_max, int mr,
                      double tol_abs, double tol_rel, int sorted=1, int blocked=0)

*/
/******************************************************************************//*


{
 double av;
 double *c;
 double delta = 1.0e-03;
 double *g;
 double **h;
 double htmp;
 int i;
 int itr;
 int itr_used;
 int j;
 int k;
 int k_copy;
 double *l;
 double mu;
 double *r;
 double rho;
 double rho_tol;
 double *s;
 int *ua;
 double **v;
 int verbose = 0;
 double *y;
 double alp[2] = {1.0, 0};
 double bet[2] = {0.0, 0};
 itr_used = 0;

 c = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
 g = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
 h = dmatrix ( 0, mr, 0, mr-1 );
 r = ( double * ) malloc ( n * sizeof ( double ) );
 s = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
 ua = ( int * ) malloc ( n * sizeof ( int ) );
 v = dmatrix ( 0, mr, 0, n-1 );
 y = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
 size_t *ap = ( size_t * ) malloc ( (n+1) * sizeof ( size_t ) );
 for (int m = 0; m < n+1; ++m) {
  ap[m] = ia[m];
  //printf("%ld == %d, ", ap[m], ia[m]);
 }
 int xtype = CHOLMOD_REAL;
 cholmod_dense *B_o = cholmod_l_zeros (n, 1, xtype, cm);
 cholmod_dense *Y ;
 double *Yx;
 double *B_ox = (double*) B_o->x;
 if (!sorted){
  rearrange_cr ( n, nz_num, ia, ja, a ); //FIXME:
  diagonal_pointer_cr ( n, nz_num, ia, ja, ua );
 }else{
  for (int m = 0; m < n; ++m) {
   ua[m] = ia[m];
  }
 }
 if ( verbose ) {
  printf ( "\n" );
  printf ( "PMGMRES_ILU_CR\n" );
  printf ( "  Number of unknowns = %d\n", n );
 }

 for ( itr = 0; itr < itr_max; itr++ ) {
  if(!sorted)
   ax_cr ( n, nz_num, ia, ja, a, x, r );
  else
   spmv_csc_sym_one(n, ap, ja, a, -1, alp, bet, 1,x, r );

*/
/*  printf("\n");
  for (int m = 0; m < n; ++m) {
   printf("%f, ",x[m]);
  }
  printf("\n");
  for (int m = 0; m < n; ++m) {
   printf("%f, ",r[m]);
  }
  printf("\n");*//*

  for ( i = 0; i < n; i++ ) {
   r[i] = rhs[i] - r[i];
   B_ox[i] =  r[i];
  }
 // print_vec("R:\n",0,n,B_ox);
*/
/*  if (blocked==2)
   solve_phase_ldl_blocked_parallel(n, d_val, r, col2sup, sup2col,
                                    Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                    levels, levelPtr, levelSet,
                                    parts, parPtr, partition, chunk);
  else if(blocked==1)
   solve_phase_ldl_blocked(n, d_val, r, col2sup, sup2col,
                           Lp, Li, Lx, Li_ptr, supNo, NNZ);
  else
   solve_phase_ldl(n, d_val, r, col2sup, sup2col,
                   Lp, Li, Lx, Li_ptr, supNo, NNZ);*//*

  Y = cholmod_l_solve(1,L,B_o,cm);
  Yx = (double *)Y->x;
  for (int m = 0; m < n; ++m) {
   r[m] = Yx[i];
  }
  cholmod_l_free_dense (&Y, cm) ;
  //print_vec("R:\n",0,n,r);
  rho = sqrt ( r8vec_dot ( n, r, r ) );
  if(rho < tol_abs)
   return itr_used;
  if ( verbose ) {
   printf ( "  ITR = %d  Residual = %e\n", itr, rho );
  }
  if ( itr == 0 ) {
   rho_tol = rho * tol_rel;
  }

  for ( i = 0; i < n; i++ ) {
   v[0][i] = r[i] / rho;
  }
  g[0] = rho;
  for ( i = 1; i < mr + 1; i++ ) {
   g[i] = 0.0;
  }
  for ( i = 0; i < mr + 1; i++ ) {
   for ( j = 0; j < mr; j++ ) {
    h[i][j] = 0.0;
   }
  }
  for ( k = 0; k < mr; k++ ) {
   k_copy = k;
   if(!sorted)
    ax_cr ( n, nz_num, ia, ja, a, v[k], v[k+1] );
   else
    spmv_csc_sym_one(n, ap, ja, a, -1, alp, bet, 1, v[k], v[k+1] );

*/
/*   if (blocked==2)
    solve_phase_ldl_blocked_parallel(n, d_val, v[k+1], col2sup, sup2col,
                                     Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                     levels, levelPtr, levelSet,
                                     parts, parPtr, partition, chunk);
   else if(blocked==1)
    solve_phase_ldl_blocked(n, d_val, v[k+1], col2sup, sup2col,
                            Lp, Li, Lx, Li_ptr, supNo, NNZ);
   else
    solve_phase_ldl(n, d_val, v[k+1], col2sup, sup2col,
                    Lp, Li, Lx, Li_ptr, supNo, NNZ);*//*

   for ( i = 0; i < n; i++ ) {
    B_ox[i] = v[k+1][i];
   }
   Y = cholmod_l_solve(1,L,B_o,cm);
   Yx = (double *)Y->x;
   for (int m = 0; m < n; ++m) {
    v[k+1][m] = Yx[i];
   }
   cholmod_l_free_dense (&Y, cm) ;

   av = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );
   for ( j = 0; j <= k; j++ ) {
    h[j][k] = r8vec_dot ( n, v[k+1], v[j] );
    for ( i = 0; i < n; i++ ) {
     v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
    }
   }
   h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

   if ( ( av + delta * h[k+1][k]) == av ) {
    for ( j = 0; j < k + 1; j++ ) {
     htmp = r8vec_dot ( n, v[k+1], v[j] );
     h[j][k] = h[j][k] + htmp;
     for ( i = 0; i < n; i++ ) {
      v[k+1][i] = v[k+1][i] - htmp * v[j][i];
     }
    }
    h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );
   }

   if ( h[k+1][k] != 0.0 ) {
    for ( i = 0; i < n; i++ ) {
     v[k+1][i] = v[k+1][i] / h[k+1][k];
    }
   }

   if ( 0 < k ) {
    for ( i = 0; i < k + 2; i++ ) {
     y[i] = h[i][k];
    }
    for ( j = 0; j < k; j++ ) {
     mult_givens ( c[j], s[j], j, y );
    }
    for ( i = 0; i < k + 2; i++ ) {
     h[i][k] = y[i];
    }
   }
   mu = sqrt ( h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k] );
   if(mu ==0)
    return itr_used;
   c[k] = h[k][k] / mu;
   s[k] = -h[k+1][k] / mu;
   h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
   h[k+1][k] = 0.0;
   mult_givens ( c[k], s[k], k, g );

   rho = fabs ( g[k+1] );

   itr_used = itr_used + 1;

   if ( verbose ) {
    printf ( "  K   = %d  Residual = %e\n", k, rho );
   }

   if ( rho <= rho_tol && rho <= tol_abs ) {
    break;
   }
  }

  k = k_copy;

  y[k] = g[k] / h[k][k];
  for ( i = k - 1; 0 <= i; i-- ) {
   y[i] = g[i];
   for ( j = i + 1; j < k + 1; j++ ) {
    y[i] = y[i] - h[i][j] * y[j];
   }
   y[i] = y[i] / h[i][i];
  }
  for ( i = 0; i < n; i++ ) {
   for ( j = 0; j < k + 1; j++ ) {
    x[i] = x[i] + v[j][i] * y[j];
   }
  }
  if ( rho <= rho_tol && rho <= tol_abs ) {
   break;
  }
 }

 if ( verbose )
 {
  printf ( "\n" );
  printf ( "PMGMRES_ILU_CR:\n" );
  printf ( "  Iterations = %d\n", itr_used );
  printf ( "  Final residual = %e\n", rho );
 }
 free ( c );
 free ( g );
 free_dmatrix ( h, 0, mr, 0, mr - 1 );
 free ( r );
 free ( s );
 free_dmatrix ( v, 0, mr, 0, n - 1 );
 free ( y );
 return itr_used;
}
*/

 void print_cs(cholmod_sparse *A) {
  long *Ap = (long *) A->p;
  long *Ai = (long *) A->i;
  double *Ax = (double *) A->x;
  for (int i = 0; i < A->ncol; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    std::cout << i << " " << Ai[j] << " " << Ax[j] << "\n";
   }
  }
  std::cout << "\n";
 }

 cholmod_sparse *expand_matrix(int k, cholmod_sparse *A, double diag,
                               cholmod_common *cm) {
  size_t n = A->ncol;
  long *Ap = (long *) (A->p);
  long *Ai = (long *) (A->i);
  double *Ax = (double *) (A->x);
  size_t new_n = n + k;
  size_t new_nnz = A->nzmax + k;
  cholmod_sparse *B = cholmod_l_allocate_sparse(new_n, new_n, new_nnz,
                                                TRUE, TRUE, 0,
                                                CHOLMOD_REAL, cm);
  long *Bp = (long *) (B->p);
  long *Bi = (long *) (B->i);
  double *Bx = (double *) (B->x);

  for (int j = 0; j < n + 1; ++j) {
   Bp[j] = Ap[j];
  }
  for (int i = n + 1; i < new_n + 1; ++i) {
   Bp[i] = Bp[i - 1] + 1;
  }

  for (int l = 0; l < A->nzmax; ++l) {
   Bi[l] = Ai[l];
   Bx[l] = Ax[l];
  }
  long cur_cul = n;
  for (int m = A->nzmax; m < new_nnz; ++m) {
   Bi[m] = cur_cul;
   Bx[m] = diag;
   cur_cul++;
  }
  return B;
 }

 void add_row_to_ordered(CSC *A, CSC *BT, size_t *iperm, int row_idx,
                         CSC *eA) {
  size_t nz_e_row = BT->p[row_idx + 1] - BT->p[row_idx];
  eA->nzmax = A->nzmax + nz_e_row;
  eA->ncol = A->ncol + 1;
  eA->nrow = A->nrow + 1;
  eA->p = new int[eA->ncol + 1];
  eA->i = new int[eA->nzmax];
  eA->x = new double[eA->nzmax];
  for (int i = 0; i < A->ncol + 1; ++i) {
   eA->p[i] = A->p[i];
  }
  eA->p[eA->ncol] = eA->p[eA->ncol - 1] + nz_e_row;
  for (int j = 0; j < A->nzmax; ++j) {
   eA->i[j] = A->i[j];
   eA->x[j] = A->x[j];
  }
  for (int k = eA->p[eA->ncol - 1], t = BT->p[row_idx]; k < nz_e_row; ++k, ++t) {
   eA->x[k] = BT->x[t];
   eA->i[k] = iperm[A->i[t]];
  }
 }

 int edit_l_factor(int row_of_kkt, int row_of_BT, cholmod_factor *L,
                   CSC *BT,
                   int upordown, cholmod_common *cm, double eps_pert) {
  size_t n = L->n;
  int st = 1;
  long *LIPerm = (long *) L->IPerm;
  long *LPerm = (long *) L->Perm;
/* for (int j = 0; j < n; ++j) {
  LIPerm[ LPerm[j] ] = j;
 }*/
  size_t kDelAdd = LIPerm[row_of_kkt];
  assert(kDelAdd < n);
  if (upordown == 0) { // downdate L
   //kDelAdd = LIPerm[row_idx];
   //int st = cholmod_l_rowdel_solve(kDelAdd,NULL,yk,L1,X,B_o,cm);
   st = cholmod_l_rowdel(kDelAdd, NULL, L, cm);
   if (!st)
    std::cout << "deletion failed\n";
  } else { // update L
   // creating the extra row from B and row_idx
   size_t base = BT->nrow;
   //kDelAdd = LIPerm[row_idx];
   size_t nnz_extra = BT->p[row_of_BT + 1] - BT->p[row_of_BT] + 1;// 1 for diagonal
   cholmod_sparse *extraRow = cholmod_l_allocate_sparse(n, 1, nnz_extra,
                                                        TRUE, TRUE, 0,
                                                        CHOLMOD_REAL, cm);
   long *erP = (long *) (extraRow->p);
   long *erI = (long *) (extraRow->i);
   double *erX = (double *) (extraRow->x);
   erP[0] = 0;
   erP[1] = nnz_extra;
   for (int i1 = BT->p[row_of_BT], i2 = 0; i1 < BT->p[row_of_BT + 1]; ++i1) {
    erI[i2] = LIPerm[BT->i[i1]];
    erX[i2] = BT->x[i1];
    i2++;
   }
   erI[nnz_extra - 1] = kDelAdd;//diagonal
   erX[nnz_extra - 1] = eps_pert;//diagonal
   // Add it to L
   st = cholmod_l_rowadd(kDelAdd, extraRow, L, cm);
   if (!st)
    std::cout << "adding failed in row: " << kDelAdd << "\n";
   cholmod_free_sparse(&extraRow, cm);
  }

  return st;
 }


 cholmod_sparse *converting_to_cs(CSC *A, cholmod_common *cm) {
  size_t n = A->ncol;
  size_t nnz = A->nzmax;
  cholmod_sparse *A_cs = cholmod_l_allocate_sparse(n, n, nnz,
                                                   TRUE, TRUE, -1, CHOLMOD_REAL, cm);
  long *A_csi = (long *) A_cs->i;
  long *A_csp = (long *) A_cs->p;
  double *A_csx = (double *) A_cs->x;
  for (int i = 0; i < n; ++i) {
   A_csp[i] = A->p[i];
  }
  A_csp[n] = A->p[n];
  for (int i = 0; i < n; ++i) {
   for (int j = A_csp[i]; j < A_csp[i + 1]; ++j) {
    A_csi[j] = A->i[j];
    A_csx[j] = A->x[j];
   }
  }
  cholmod_sparse *C = cholmod_l_transpose(A_cs, 2, cm);
  cholmod_l_free_sparse(&A_cs, cm);
  return C;
 }


 void chol_solve_updated_kkt(double *rhs, cholmod_factor *L,
                             cholmod_common *cm, double *x, CSC *kkt_updated, CSC *H, CSC *B, CSC *C,
                             std::vector<int> as,
                             int outer_iter) {
  size_t n = L->n;
  int xtype = L->xtype;
  long *LPerm = (long *) L->Perm;
  cholmod_dense *B_o = cholmod_l_zeros(n, 1, xtype, cm);
  double *B_ox = (double *) B_o->x;
  for (int k1 = 0; k1 < n; ++k1) {
   B_ox[k1] = rhs[LPerm[k1]];
  }
  //print_vec("B:\n",0,n,B_ox);
  cholmod_dense *Y;
  Y = cholmod_l_solve(1, L, B_o, cm);
  double *Yx = (double *) Y->x;
  //print_vec("X:\n",0,n,Yx);

/* cholmod_dense *YY, *ZZ = cholmod_l_ones(n,1,CHOLMOD_REAL,cm) ;
 YY = cholmod_l_solve(1,L,ZZ,cm);
 double *YYx = (double*)YY->x;
 print_vec("X:\n",0,n,YYx);*/
  /// Iterative refinement
  if (outer_iter > 0) {
   int *perm_int = new int[n]();
   for (int j = 0; j < n; ++j) {
    perm_int[j] = LPerm[j];
   }

   CSC *A_red = new CSC;
   CSC *kkt_u_t = new CSC;
   CSC *kkt_u_t_tmp = new CSC;
   double *rhs_test;
   gather_rows(as, B->nrow, B->ncol, B->p, B->i, B->x,
               A_red->nrow, A_red->ncol, A_red->nzmax,
               A_red->p, A_red->i, A_red->x);
   build_super_matrix_with_eq_chol(n, H, A_red, C, 0, kkt_u_t->ncol,
                                   kkt_u_t->nzmax,
                                   kkt_u_t->p, kkt_u_t->i, kkt_u_t->x,
                                   rhs_test);
   kkt_u_t->nrow = kkt_u_t->ncol;
   kkt_updated->p = new int[kkt_u_t->ncol + 1]();
   kkt_updated->i = new int[kkt_u_t->nzmax]();
   kkt_updated->x = new double[kkt_u_t->nzmax]();
   kkt_updated->ncol = kkt_updated->nrow = kkt_u_t->ncol;
   kkt_updated->nzmax = kkt_u_t->p[kkt_u_t->ncol];
   kkt_updated->stype = -1;
   kkt_updated->packed = 1;
   kkt_updated->xtype = CHOLMOD_REAL;
   kkt_updated->sorted = 1;

   kkt_u_t_tmp->p = new int[kkt_u_t->ncol + 1]();
   kkt_u_t_tmp->i = new int[kkt_u_t->nzmax]();
   kkt_u_t_tmp->x = new double[kkt_u_t->nzmax]();
   kkt_u_t_tmp->ncol = kkt_u_t_tmp->nrow = kkt_u_t->ncol;
   kkt_u_t_tmp->nzmax = kkt_u_t->p[kkt_u_t->ncol];
   kkt_u_t_tmp->stype = 1;
   kkt_u_t_tmp->packed = 1;
   kkt_u_t_tmp->xtype = CHOLMOD_REAL;
   kkt_u_t_tmp->sorted = 1;

   kkt_u_t->stype = -1;
   kkt_u_t->packed = 1;
   kkt_u_t->xtype = CHOLMOD_REAL;
   kkt_u_t->sorted = 1;
   int ret_val = transpose_sym(kkt_u_t, 1, perm_int, kkt_u_t_tmp, 1);
   transpose_sym(kkt_u_t_tmp, 1, perm_int, kkt_updated, 1);
   assert(kkt_updated->ncol == n);
   pmgmres_ldlt_auto_chol(kkt_updated->ncol, kkt_updated->nzmax, kkt_updated->p,
                          kkt_updated->i,
                          kkt_updated->x, L, cm, Yx, B_ox, outer_iter, outer_iter, 1e-15, 1e-15);
   allocateAC(A_red, 0, 0, 0, FALSE);
   allocateAC(kkt_u_t, 0, 0, 0, FALSE);
   allocateAC(kkt_u_t_tmp, 0, 0, 0, FALSE);
   allocateAC(kkt_updated, 0, 0, 0, FALSE);
   delete[]rhs_test;
  }

  for (int i = 0; i < n; ++i) {
   x[LPerm[i]] = Yx[i];
  }
  cholmod_free_dense(&B_o, cm);
  cholmod_free_dense(&Y, cm);

 }

 cholmod_factor *initial_solve_1(cholmod_sparse *kkt, CSC *A,
                                 CSC *C, //eq
                                 int num_ineq, double *rhs, double *x,
                                 cholmod_common *cm, double eps) {
// Creating KKT from hessian, equalities and add enough diagonals
  size_t nnz_k = A->nzmax + (C ? C->nzmax : 0) + (C ? C->nrow : 0) + num_ineq;
  size_t n = A->ncol + num_ineq + (C ? C->nrow : 0);
  kkt = cholmod_l_allocate_sparse(n, n, nnz_k, TRUE, TRUE, -1,
                                  CHOLMOD_REAL, cm);
  cholmod_factor *L;
  long *KP = (long *) (kkt->p);
  long *KI = (long *) (kkt->i);
  double *KX = (double *) (kkt->x);
  KP[0] = 0;
  if (C ? C->nrow > 0 : 0) {
   for (int i = 1; i < A->ncol + 1; ++i) {
    KP[i] = KP[i - 1] + (A->p[i] - A->p[i - 1]) +
            (C->p[i] - C->p[i - 1]);
   }
  } else {
   for (int i = 1; i < A->ncol + 1; ++i) {
    KP[i] = KP[i - 1] + (A->p[i] - A->p[i - 1]);
   }
  }
  //Adding diagonal for columns with zero values, both eq and ineq.
  for (int k = A->ncol + 1; k < n + 1; ++k) {
   KP[k] = KP[k - 1] + 1;
  }
  assert(KP[n] == nnz_k);
  int base1 = A->ncol;
  size_t stp = 0;
  for (int j = 0; j < A->ncol; ++j) {
   stp = KP[j];
   //Adding Hessian
   for (int i = A->p[j]; i < A->p[j + 1]; ++i) {
    KI[stp] = A->i[i];
    KX[stp] = A->x[i];
    stp++;
   }
   base1 = A->ncol;
   //Adding equalities
   if (C != NULL) {
    if (C->nrow > 0) {
     for (int i = C->p[j]; i < C->p[j + 1]; ++i) {
      KI[stp] = base1 + C->i[i];
      KX[stp] = C->x[i];
      //std::cout<<"Eq: "<< base1 + C->i[i]<<"; "<<SMx[stp]<<"\n";
      stp++;
     }
    }
   }
  }
  size_t base2 = A->ncol;
  for (int l = base2; l < n; ++l) {
   for (int i = KP[l]; i < KP[l + 1]; ++i) {
    KI[i] = base2;
    KX[i] = eps;
    base2++;
   }
  }
  // factorize and solve
  cholmod_dense *X, *BB;
  L = cholmod_l_analyze(kkt, cm);
  int xtype = CHOLMOD_REAL;
  cholmod_dense *B_rhs = cholmod_l_zeros(n, 1, xtype, cm);
  double *B_rhsx = (double *) (B_rhs->x);
  for (int k = 0; k < n; ++k) {
   B_rhsx[k] = rhs[k];
  }
  //print_cs(kkt);
  //print_vec("\nrhs:\n",0,n,B_rhsx);
  cholmod_l_factorize(kkt, L, cm);
  X = cholmod_l_solve(CHOLMOD_A, L, B_rhs, cm);
  double *Xx = (double *) X->x;
// copy solution to x
  for (int l = 0; l < n; ++l) {
   x[l] = Xx[l];
  }
  long *LPerm = (long *) L->Perm;
  L->IPerm = (long *) calloc(n, sizeof(long));
  long *LIPerm = (long *) L->IPerm;
  for (int j = 0; j < n; ++j) {
   LIPerm[LPerm[j]] = j;
  }
  // reorder the kkt system for iterative refinement
/* base2 = A->ncol ;
 //reset diagonals
 for (int l = base2; l < n; ++l) {
  for (int i = KP[l]; i < KP[l + 1]; ++i) {
   KX[i] = 0;
  }
 }
 cholmod_sparse *kkt_tmp=NULL;
 kkt_tmp = cholmod_l_ptranspose(kkt,1,LPerm, NULL, 0, cm);
 kktT = cholmod_l_ptranspose(kkt_tmp,1,NULL, NULL, 0, cm);
 cholmod_l_free_sparse(&kkt_tmp,cm);

 for (int l = base2; l < n; ++l) {
  for (int i = KP[l]; i < KP[l + 1]; ++i) {
   KX[i] = eps;
  }
 }*/
  return L;
 }

 void initial_solve_2(cholmod_sparse *kkt, cholmod_factor *L, CSC *inc_mat,
                      double *rhs, double *x, cholmod_common *cm) {
  size_t n = inc_mat->ncol;
  size_t nnz_k = inc_mat->nzmax;
  int xtype = L->xtype;
  cholmod_dense *X;
// Copying super KKT to cholmod_sparse
  kkt = cholmod_l_allocate_sparse(n, 1, nnz_k, TRUE, TRUE, 0, CHOLMOD_REAL, cm);
  long *KP = (long *) (kkt->p);
  long *KI = (long *) (kkt->i);
  double *KX = (double *) (kkt->x);
  for (int i = 0; i < n; ++i) {
   KP[i] = inc_mat->p[i];
   for (int j = inc_mat->p[i]; j < inc_mat->p[i + 1]; ++j) {
    KI[j] = inc_mat->i[j];
    KX[j] = inc_mat->x[j];
   }
  }
  KP[n] = inc_mat->p[n];
  cholmod_dense *B = cholmod_l_zeros(n, 1, xtype, cm);
  double *Bx = (double *) (B->x);
  for (int k = 0; k < n; ++k) {
   Bx[k] = rhs[k];
  }
// factorize and solve
  L = cholmod_l_analyze(kkt, cm);
  cholmod_l_factorize(kkt, L, cm);
  X = cholmod_l_solve(CHOLMOD_A, L, B, cm);
  double *Xx = (double *) X->x;
// copy solution to x
  for (int l = 0; l < n; ++l) {
   x[l] = Xx[l];
  }
 }


}

#endif //PARS_CHOLMOD_UTILS_H
