//
// Created by Shujian Qian on 2020-10-29.
//

#include "nasoq/gmres/mgmres.hpp"

# include <cstdio>
# include <cstdlib>
# include <ctime>

#include "nasoq/common/Norm.h"
#include "nasoq/linear_solver/solve_phase.h"
#include "nasoq/matrixVector/spmv_CSC.h"

namespace nasoq {
 /******************************************************************************/

 void ax_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[],
            double w[])

/******************************************************************************/
/*
  Purpose:

    AX_CR computes A*x for a matrix stored in sparse compressed row form.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    For this version of MGMRES, the row and column indices are assumed
    to use the C/C++ convention, in which indexing begins at 0.

    If your index vectors IA and JA are set up so that indexing is based
    at 1, each use of those vectors should be shifted down by 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double A[NZ_NUM], the matrix values.

    Input, double X[N], the vector to be multiplied by A.

    Output, double W[N], the value of A*X.
*/
 {
  int i;
  int k;
  int k1;
  int k2;

  for (i = 0; i < n; i++) {
   w[i] = 0.0;
   k1 = ia[i];
   k2 = ia[i + 1];
   for (k = k1; k < k2; k++) {
    w[i] = w[i] + a[k] * x[ja[k]];
   }
  }
  return;
 }
/******************************************************************************/


/******************************************************************************/

 void diagonal_pointer_cr(int n, int nz_num, int ia[], int ja[], int ua[])

/******************************************************************************/
/*
  Purpose:

    DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    The array UA can be used to locate the diagonal elements of the matrix.

    It is assumed that every row of the matrix includes a diagonal element,
    and that the elements of each row have been ascending sorted.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.  On output,
    the order of the entries of JA may have changed because of the sorting.

    Output, int UA[N], the index of the diagonal element of each row.
*/
 {
  int i;
  int j;
  int j1;
  int j2;
  int k;

  for (i = 0; i < n; i++) {
   ua[i] = -1;
   j1 = ia[i];
   j2 = ia[i + 1];

   for (j = j1; j < j2; j++) {
    if (ja[j] == i) {
     ua[i] = j;
    }
   }

  }
  return;
 }

 /******************************************************************************/

 void mult_givens(double c, double s, int k, double *g)

/******************************************************************************/
/*
  Purpose:

    MULT_GIVENS applies a Givens rotation to two vector elements.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, double C, S, the cosine and sine of a Givens
    rotation.

    Input, int K, indicates the location of the first vector entry.

    Input/output, double G[K+2], the vector to be modified.  On output,
    the Givens rotation has been applied to entries G(K) and G(K+1).
*/
 {
  double g1;
  double g2;

  g1 = c * g[k] - s * g[k + 1];
  g2 = s * g[k] + c * g[k + 1];

  g[k] = g1;
  g[k + 1] = g2;

  return;
 }

 /******************************************************************************/

 double r8vec_dot(int n, double a1[], double a2[])

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT, the dot product of the vectors.
*/
 {
  int i;
  double value;

  value = 0.0;
  for (i = 0; i < n; i++) {
   value = value + a1[i] * a2[i];
  }
  return value;
 }

 /******************************************************************************/

 void rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[])

/******************************************************************************/
/*
  Purpose:

    REARRANGE_CR sorts a sparse compressed row matrix.

  Discussion:

    This routine guarantees that the entries in the CR matrix
    are properly sorted.

    After the sorting, the entries of the matrix are rearranged in such
    a way that the entries of each column are listed in ascending order
    of their column values.

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], the compressed row index vector.

    Input/output, int JA[NZ_NUM], the column indices of the matrix values.
    On output, the order of the entries of JA may have changed because of
    the sorting.

    Input/output, double A[NZ_NUM], the matrix values.  On output, the
    order of the entries may have changed because of the sorting.
*/
 {
  double dtemp;
  int i;
  int is;
  int itemp;
  int j;
  int j1;
  int j2;
  int k;

  for (i = 0; i < n; i++) {
   j1 = ia[i];
   j2 = ia[i + 1];
   is = j2 - j1;

   for (k = 1; k < is; k++) {
    for (j = j1; j < j2 - k; j++) {
     if (ja[j + 1] < ja[j]) {
      itemp = ja[j + 1];
      ja[j + 1] = ja[j];
      ja[j] = itemp;

      dtemp = a[j + 1];
      a[j + 1] = a[j];
      a[j] = dtemp;
     }
    }
   }
  }
  return;
 }

/******************************************************************************/

 void timestamp()

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
 {
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
# undef TIME_SIZE
 }

 double **dmatrix(int nrl, int nrh, int ncl, int nch)

/******************************************************************************/
/*
  Purpose:

    DMATRIX allocates a double matrix.

  Discussion:

    The matrix will have a subscript range m[nrl...nrh][ncl...nch] .

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Output, double **DMATRIX, a doubly-dimensioned array with
    the requested row and column ranges.
*/
 {
  int i;
  double **m;
  int nrow = nrh - nrl + 1;
  int ncol = nch - ncl + 1;
/*
  Allocate pointers to the rows.
*/
  m = (double **) malloc((size_t) ((nrow + 1) * sizeof(double *)));

  if (!m) {
   fprintf(stderr, "\n");
   fprintf(stderr, "DMATRIX - Fatal error!\n");
   fprintf(stderr, "  Failure allocating pointers to rows.\n");
   exit(1);
  }
  m = m + 1;
  m = m - nrl;
/*
  Allocate each row and set pointers to them.
*/
  m[nrl] = (double *) malloc((size_t) ((nrow * ncol + 1) * sizeof(double)));

  if (!m[nrl]) {
   fprintf(stderr, "\n");
   fprintf(stderr, "DMATRIX - Fatal error!\n");
   fprintf(stderr, "  Failure allocating rows.\n");
   exit(1);
  }
  m[nrl] = m[nrl] + 1;
  m[nrl] = m[nrl] - ncl;

  for (i = nrl + 1; i <= nrh; i++) {
   m[i] = m[i - 1] + ncol;
  }
/*
  Return the pointer to the array of pointers to the rows;
*/
  return m;
 }

 double **dmatrix_prealloc(int nrl, int nrh, int ncl, int nch, double *ws)

/******************************************************************************/
/*

*/
 {
  int i;
  double **m;
  int nrow = nrh - nrl + 1;
  int ncol = nch - ncl + 1;
/*
  Allocate pointers to the rows.
*/
  m = (double **) malloc((size_t) ((nrow + 1) * sizeof(double *)));

  if (!m) {
   fprintf(stderr, "\n");
   fprintf(stderr, "DMATRIX - Fatal error!\n");
   fprintf(stderr, "  Failure allocating pointers to rows.\n");
   exit(1);
  }
  m = m + 1;
  m = m - nrl;
/*
  Allocate each row and set pointers to them.
*/
  //m[nrl] = ( double * ) malloc ( (size_t) ( ( nrow * ncol + 1 ) * sizeof ( double ) ) );
  m[nrl] = ws;
  if (!m[nrl]) {
   fprintf(stderr, "\n");
   fprintf(stderr, "DMATRIX - Fatal error!\n");
   fprintf(stderr, "  Failure allocating rows.\n");
   exit(1);
  }
  m[nrl] = m[nrl] + 1;
  m[nrl] = m[nrl] - ncl;

  for (i = nrl + 1; i <= nrh; i++) {
   m[i] = m[i - 1] + ncol;
  }
/*
  Return the pointer to the array of pointers to the rows;
*/
  return m;
 }

 void free_dmatrix_prealloc(double **m, int nrl, int nrh, int ncl, int nch)

/******************************************************************************/
/*

*/
 {
  //free ( ( char * ) ( m[nrl] + ncl - 1 ) );
  free((char *) (m + nrl - 1));
  //free(m);
  return;
 }

 void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)

/******************************************************************************/
/*
  Purpose:

    FREE_DMATRIX frees a double matrix allocated by DMATRIX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Input, double **M, the pointer to the doubly-dimensioned array,
    previously created by a call to DMATRIX.
*/
 {
  free((char *) (m[nrl] + ncl - 1));
  free((char *) (m + nrl - 1));

  return;
 }

 void lus_cr(int n, int nz_num, int *ia, int *ja, double *l, int *ua, double *r, double *z)

/******************************************************************************/
/*
  Purpose:

    LUS_CR applies the incomplete LU preconditioner.

  Discussion:

    The linear system M * Z = R is solved for Z.  M is the incomplete
    LU preconditioner matrix, and R is a vector supplied by the user.
    So essentially, we're solving L * U * Z = R.

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double L[NZ_NUM], the matrix values.

    Input, int UA[N], the index of the diagonal element of each row.

    Input, double R[N], the right hand side.

    Output, double Z[N], the solution of the system M * Z = R.
*/
 {
  int i;
  int j;
  double *w;

  w = (double *) malloc(n * sizeof(double));
/*
  Copy R in.
*/
  for (i = 0; i < n; i++) {
   w[i] = r[i];
  }
/*
  Solve L * w = w where L is unit lower triangular.
*/
  for (i = 1; i < n; i++) {
   for (j = ia[i]; j < ua[i]; j++) {
    w[i] = w[i] - l[j] * w[ja[j]];
   }
  }
/*
  Solve U * w = w, where U is upper triangular.
*/
  for (i = n - 1; 0 <= i; i--) {
   for (j = ua[i] + 1; j < ia[i + 1]; j++) {
    w[i] = w[i] - l[j] * w[ja[j]];
   }
   w[i] = w[i] / l[ua[i]];
  }
/*
  Copy Z out.
*/
  for (i = 0; i < n; i++) {
   z[i] = w[i];
  }
/*
  Free memory.
*/
  free(w);

  return;
 }

 int pmgmres_ldlt_cr(int n, int nz_num, int *ia, int *ja, double *a, size_t *Lp, int *Li, double *Lx, size_t NNZ,
                     size_t *Li_ptr, int *col2sup, int *sup2col, int supNo, double *d_val, double *x, double *rhs,
                     int itr_max, int mr, double tol_abs, double tol_rel, int sorted, int blocked, int levels,
                     int *levelPtr, int *levelSet, int parts, int *parPtr, int *partition, int chunk)

/******************************************************************************/
/*
  Purpose:

    PMGMRES_IDL_CR applies the preconditioned restarted GMRES algorithm.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    This routine uses the incomplete LU decomposition for the
    preconditioning.  This preconditioner requires that the sparse
    matrix data structure supplies a storage position for each diagonal
    element of the matrix A, and that each diagonal element of the
    matrix A is not zero.

    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
    corrections to the code on 31 May 2007.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994.
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 2003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the linear system.

    Input, int NZ_NUM, the number of nonzero matrix values.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double A[NZ_NUM], the matrix values.

    Input/output, double X[N]; on input, an approximation to
    the solution.  On output, an improved approximation.

    Input, double RHS[N], the right hand side of the linear system.

    Input, int ITR_MAX, the maximum number of (outer) iterations to take.

    Input, int MR, the maximum number of (inner) iterations to take.
    MR must be less than N.

    Input, double TOL_ABS, an absolute tolerance applied to the
    current residual.

    Input, double TOL_REL, a relative tolerance comparing the
    current residual to the initial residual.
*/
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

  c = (double *) malloc((mr + 1) * sizeof(double));
  g = (double *) malloc((mr + 1) * sizeof(double));
  h = dmatrix(0, mr, 0, mr - 1);
  l = (double *) malloc((ia[n] + 1) * sizeof(double));
  r = (double *) malloc(n * sizeof(double));
  s = (double *) malloc((mr + 1) * sizeof(double));
  ua = (int *) malloc(n * sizeof(int));
  v = dmatrix(0, mr, 0, n - 1);
  y = (double *) malloc((mr + 1) * sizeof(double));
  size_t *ap = (size_t *) malloc((n + 1) * sizeof(size_t));
  for (int m = 0; m < n + 1; ++m) {
   ap[m] = ia[m];
   //printf("%ld == %d, ", ap[m], ia[m]);
  }

  if (!sorted) {
   rearrange_cr(n, nz_num, ia, ja, a); //FIXME:
   diagonal_pointer_cr(n, nz_num, ia, ja, ua);
  } else {
   for (int m = 0; m < n; ++m) {
    ua[m] = ia[m];
   }
  }
  if (verbose) {
   printf("\n");
   printf("PMGMRES_ILU_CR\n");
   printf("  Number of unknowns = %d\n", n);
  }

  for (itr = 0; itr < itr_max; itr++) {
   if (!sorted)
    ax_cr(n, nz_num, ia, ja, a, x, r);
   else
    spmv_csc_sym_one(n, ap, ja, a, -1, alp, bet, 1, x, r);

/*  printf("\n");
  for (int m = 0; m < n; ++m) {
   printf("%f, ",x[m]);
  }
  printf("\n");
  for (int m = 0; m < n; ++m) {
   printf("%f, ",r[m]);
  }
  printf("\n");*/
   for (i = 0; i < n; i++) {
    r[i] = rhs[i] - r[i];
   }
   if (blocked == 2)
    solve_phase_ldl_blocked_parallel(n, d_val, r, col2sup, sup2col,
                                     Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                     levels, levelPtr, levelSet,
                                     parts, parPtr, partition, chunk);
   else if (blocked == 1)
    solve_phase_ldl_blocked(n, d_val, r, col2sup, sup2col,
                            Lp, Li, Lx, Li_ptr, supNo, NNZ);
   else
    solve_phase_ldl(n, d_val, r, col2sup, sup2col,
                    Lp, Li, Lx, Li_ptr, supNo, NNZ);

   rho = sqrt(r8vec_dot(n, r, r));
   if (rho < tol_abs)
    return itr_used;
   if (verbose) {
    printf("  ITR = %d  Residual = %e\n", itr, rho);
   }
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
     spmv_csc_sym_one(n, ap, ja, a, -1, alp, bet, 1, v[k], v[k + 1]);

    if (blocked == 2)
     solve_phase_ldl_blocked_parallel(n, d_val, v[k + 1], col2sup, sup2col,
                                      Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                      levels, levelPtr, levelSet,
                                      parts, parPtr, partition, chunk);
    else if (blocked == 1)
     solve_phase_ldl_blocked(n, d_val, v[k + 1], col2sup, sup2col,
                             Lp, Li, Lx, Li_ptr, supNo, NNZ);
    else
     solve_phase_ldl(n, d_val, v[k + 1], col2sup, sup2col,
                     Lp, Li, Lx, Li_ptr, supNo, NNZ);

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

    if (verbose) {
     printf("  K   = %d  Residual = %e\n", k, rho);
    }

    if (rho <= rho_tol && rho <= tol_abs) {
     break;
    }
   }

   k = k_copy;

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

  if (verbose) {
   printf("\n");
   printf("PMGMRES_ILU_CR:\n");
   printf("  Iterations = %d\n", itr_used);
   printf("  Final residual = %e\n", rho);
  }
  return itr_used;
 }

 void
 free_ws_ir(double *ws, double *c, double *g, double **h, int mr, double *r, double *s, double **v, int n, double *y) {

  if (ws == NULL) {
   free(c);
   free(g);
   free_dmatrix(h, 0, mr, 0, mr - 1);
   free(r);
   free(s);
   free_dmatrix(v, 0, mr, 0, n - 1);
   free(y);
  } else {
   //free_dmatrix ( h, 0, mr, 0, mr - 1 );
   //free_dmatrix ( v, 0, mr, 0, n - 1 );
   free_dmatrix_prealloc(h, 0, mr, 0, mr - 1);
   free_dmatrix_prealloc(v, 0, mr, 0, n - 1);
  }
 }

 int pmgmres_ldlt_auto(int n, int nz_num, int *ia, int *ja, double *a, size_t *Lp, int *Li, double *Lx, size_t NNZ,
                       size_t *Li_ptr, int *col2sup, int *sup2col, int supNo, double *d_val, double *x, double *rhs,
                       int itr_max, int mr, double tol_abs, double tol_rel, int sorted, int blocked, int levels,
                       int *levelPtr, int *levelSet, int parts, int *parPtr, int *partition, int chunk, double *ws)

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

  for (itr = 0; itr < itr_max; itr++) {
   if (!sorted)
    ax_cr(n, nz_num, ia, ja, a, x, r);
   else
    spmv_csc_sym_one_int(n, ia, ja, a, -1, alp, bet, 1, x, r);

   for (i = 0; i < n; i++) {
    r[i] = rhs[i] - r[i];
   }
   //norm_inf = norm_dense();
   if (blocked == 2)
    solve_phase_ldl_blocked_parallel(n, d_val, r, col2sup, sup2col,
                                     Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                     levels, levelPtr, levelSet,
                                     parts, parPtr, partition, chunk);
   else if (blocked == 1)
    solve_phase_ldl_blocked(n, d_val, r, col2sup, sup2col,
                            Lp, Li, Lx, Li_ptr, supNo, NNZ);
   else
    solve_phase_ldl(n, d_val, r, col2sup, sup2col,
                    Lp, Li, Lx, Li_ptr, supNo, NNZ);

   rho = sqrt(r8vec_dot(n, r, r));
   if (rho < tol_abs) {
    free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
    return itr_used;
   }
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

    if (blocked == 2)
     solve_phase_ldl_blocked_parallel(n, d_val, v[k + 1], col2sup, sup2col,
                                      Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                      levels, levelPtr, levelSet,
                                      parts, parPtr, partition, chunk);
    else if (blocked == 1)
     solve_phase_ldl_blocked(n, d_val, v[k + 1], col2sup, sup2col,
                             Lp, Li, Lx, Li_ptr, supNo, NNZ);
    else
     solve_phase_ldl(n, d_val, v[k + 1], col2sup, sup2col,
                     Lp, Li, Lx, Li_ptr, supNo, NNZ);

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
    if (mu == 0) {
     free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
     return itr_used;
    }
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
   if (h[k][k] == 0) {
    free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
    return itr_used;
   }
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
  free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
/* free ( c );
 free ( g );
 free_dmatrix ( h, 0, mr, 0, mr - 1 );
 free ( r );
 free ( s );
 free_dmatrix ( v, 0, mr, 0, n - 1 );
 free ( y );*/
  return itr_used;
 }

 int
 pmgmres_ldlt_auto_update(int n, int nz_num, int *ia, int *ja, double *a, size_t *Lp, int *Li, double *Lx, size_t NNZ,
                          size_t *Li_ptr, int *col2sup, int *sup2col, int supNo, double *d_val, double *x, double *rhs,
                          int itr_max, int mr, double tol_abs, double tol_rel, bool *mask, int *mask_col, int sorted,
                          int blocked, int levels, int *levelPtr, int *levelSet, int parts, int *parPtr, int *partition,
                          int s_level_no, int *s_level_ptr, int *s_level_set, int chunk, double *ws, double *ws_zerod)

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

  if (ws == NULL) {
   c = (double *) malloc((mr + 1) * sizeof(double));
   g = (double *) malloc((mr + 1) * sizeof(double));
   h = dmatrix(0, mr, 0, mr - 1);
   r = (double *) malloc(n * sizeof(double));
   s = (double *) malloc((mr + 1) * sizeof(double));
   v = dmatrix(0, mr, 0, n - 1);
   y = (double *) malloc((mr + 1) * sizeof(double));
  } else {
   //std::fill_n(ws,4*(mr+1) + n,0);
   // ws size is : 4*(mr+1) + mr*(mr+1) + n + mr*(n-1)
   c = ws;
   g = ws + mr + 1;
   r = ws + 2 * (mr + 1);
   s = ws + 2 * (mr + 1) + n;
   y = ws + 3 * (mr + 1) + n;
   h = dmatrix_prealloc(0, mr, 0, mr - 1, ws + 4 * (mr + 1) + n);
   //h = dmatrix ( 0, mr, 0, mr-1 );
   v = dmatrix_prealloc(0, mr, 0, n - 1, ws + 4 * (mr + 1) + mr * (mr + 1) + n);
   //v = dmatrix ( 0, mr, 0, n-1 );
  }

  for (itr = 0; itr < itr_max; itr++) {
   if (!sorted)
    ax_cr(n, nz_num, ia, ja, a, x, r);
   else
    spmv_csc_sym_one_int(n, ia, ja, a, -1, alp, bet, 1, x, r);

   for (i = 0; i < n; i++) {
    r[i] = rhs[i] - r[i];
   }
   //norm_inf = norm_dense();
   if (blocked == 2) {
    solve_phase_ldl_blocked_parallel_permuted_update(n, d_val, r, col2sup, sup2col,
                                                     Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                                     levels, levelPtr, levelSet,
                                                     parts, parPtr, partition, chunk,
                                                     s_level_no, s_level_ptr, s_level_set,
                                                     mask, mask_col, ws_zerod);
   } else if (blocked == 1)
    solve_phase_ldl_blocked(n, d_val, r, col2sup, sup2col,
                            Lp, Li, Lx, Li_ptr, supNo, NNZ);
   else
    solve_phase_ldl(n, d_val, r, col2sup, sup2col,
                    Lp, Li, Lx, Li_ptr, supNo, NNZ);

   rho = sqrt(r8vec_dot(n, r, r));
   if (rho < tol_abs) {
    free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
    return itr_used;
   }
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

    if (blocked == 2)
     solve_phase_ldl_blocked_parallel_permuted_update(n, d_val, v[k + 1], col2sup, sup2col,
                                                      Lp, Li, Lx, Li_ptr, supNo, NNZ,
                                                      levels, levelPtr, levelSet,
                                                      parts, parPtr, partition, chunk,
                                                      s_level_no, s_level_ptr, s_level_set,
                                                      mask, mask_col, ws_zerod);
    else if (blocked == 1)
     solve_phase_ldl_blocked(n, d_val, v[k + 1], col2sup, sup2col,
                             Lp, Li, Lx, Li_ptr, supNo, NNZ);
    else
     solve_phase_ldl(n, d_val, v[k + 1], col2sup, sup2col,
                     Lp, Li, Lx, Li_ptr, supNo, NNZ);

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
    if (mu == 0) {
     free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
     return itr_used;
    }
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
   if (h[k][k] == 0) {
    free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
    return itr_used;
   }
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
  free_ws_ir(ws, c, g, h, mr, r, s, v, n, y);
  return itr_used;
 }
}