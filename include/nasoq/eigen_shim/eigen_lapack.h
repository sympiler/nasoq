#pragma once

#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>
#if !EIGEN_VERSION_AT_LEAST(3, 3, 0)
#error NASOQ Eigen-based LAPACK functions require Eigen version at least 3.3.0
#endif


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
inline int LAPACKE_dsytrf(
	int matrix_layout,
	char uplo,
	int n,
	double* a,
	int lda,
	int* ipiv 
) {
	using Stride = Eigen::Stride<Eigen::Dynamic, 1>;

	if('L' != uplo && 'U' != uplo) {
		return -2; // argument 2 has an illegal value
	}

	auto copy_transpositions = [n,ipiv](const auto& ldlt) {
		const auto& P = ldlt.transpositionsP();
		assert(P.size() == n);
		for(int i = 0; i < n; ++i) ipiv[i] = static_cast<int>(1+P[i]); // `1+` converts from 0-based to 1-based
	};

	auto solve_ldlt = [uplo,n,ipiv,&copy_transpositions](auto* A_type, auto& A) {
		using Matrix = std::decay_t<decltype(*A_type)>;
		if('L' == uplo) {
			Eigen::LDLT<Eigen::Ref<Matrix>,Eigen::Lower> ldlt(A);
			if(Eigen::Success != ldlt.info()) return 1; // best we can do based on `ldlt.info()`
			copy_transpositions(ldlt);
			return 0;
		} else {
			Eigen::LDLT<Eigen::Ref<Matrix>,Eigen::Upper> ldlt(A);
			if(Eigen::Success != ldlt.info()) return 1; // best we can do based on `ldlt.info()`
			copy_transpositions(ldlt);
			return 0;
		}
	};

	if(LAPACK_ROW_MAJOR == matrix_layout) {
		using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
		Eigen::Map<Matrix, Eigen::Unaligned, Stride> A(a, n, n, Stride(lda, 1));
		return solve_ldlt((Matrix*)nullptr, A);
	} else if(LAPACK_COL_MAJOR == matrix_layout) {
		using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
		Eigen::Map<Matrix, Eigen::Unaligned, Stride> A(a, n, n, Stride(lda, 1));
		return solve_ldlt((Matrix*)nullptr, A);
	} else {
		return -1; // argument 1 has an illegal value
	}
}


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
inline int LAPACKE_dlapmt(
	int matrix_layout,
	int forwrd, // lapack_logical
	int m,
	int n,
	double* x,
	int ldx,
	int* k
) {
	using Stride = Eigen::Stride<Eigen::Dynamic, 1>;
	using PermInds = Eigen::Map<Eigen::Array<int,Eigen::Dynamic,1>>;

	auto do_permutation = [forwrd](auto& X, const auto& P) {
		if(forwrd) {
			X = X*P;
		} else {
			X = X*P.transpose();
		}
		return 0;
	};

	PermInds K_inds(k, n);
	K_inds -= 1; // 1-based to 0-based
	Eigen::PermutationWrapper<PermInds> P(K_inds);

	if(LAPACK_ROW_MAJOR == matrix_layout) {
		using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
		Eigen::Map<Matrix, Eigen::Unaligned, Stride> X(x, m, n, Stride(ldx, 1));
		do_permutation(X, P);
	} else if(LAPACK_COL_MAJOR == matrix_layout) {
		using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
		Eigen::Map<Matrix, Eigen::Unaligned, Stride> X(x, m, n, Stride(ldx, 1));
		do_permutation(X, P);
	} else {
		return -1; // argument 1 has an illegal value
	}

	K_inds += 1; // 0-based to 1-based
	return 0;
}
