#pragma once


#ifndef LAPACK_ROW_MAJOR
#define LAPACK_ROW_MAJOR 101
#endif

#ifndef LAPACK_COL_MAJOR
#define LAPACK_COL_MAJOR 102
#endif


/**
 * DSYTRF computes the factorization of a real symmetric matrix A using
 * the Bunch-Kaufman diagonal pivoting method.  The form of the
 * factorization is
 *
 *    A = U**T*D*U  or  A = L*D*L**T
 *
 * where U (or L) is a product of permutation and unit upper (lower)
 * triangular matrices, and D is symmetric and block diagonal with
 * 1-by-1 and 2-by-2 diagonal blocks.
 *
 * This is the blocked version of the algorithm, calling Level 3 BLAS.
 */
int LAPACKE_dsytrf(
	int matrix_layout,
	char uplo,
	int n,
	double* a,
	int lda,
	int* ipiv 
);


/**
 * DLAPMT rearranges the columns of the M by N matrix X as specified
 * by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
 * If FORWRD = .TRUE.,  forward permutation:
 *
 *      X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
 *
 * If FORWRD = .FALSE., backward permutation:
 *
 *      X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.
 */
int LAPACKE_dlapmt(
	int matrix_layout,
	int forwrd, // lapack_logical
	int m,
	int n,
	double* x,
	int ldx,
	int* k
);
