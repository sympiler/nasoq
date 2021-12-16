#include "nasoq/clapacke/clapacke.h"

#include <vector>

extern "C" {
#include "f2c.h"
#include "clapack.h"
}


namespace nasoq {
namespace clapacke {

template<typename T>
inline const T& minimum(const T& a, const T& b) {
	return a <= b ? a : b;
}

template<typename T>
inline const T& maximum(const T& a, const T& b) {
	return a >= b ? a : b;
}


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
	clapack_int n,
	double* a,
	clapack_int lda,
	clapack_int* ipiv 
) {
	if('L' != uplo && 'U' != uplo) {
		return -2; // argument 2 has an illegal value
	}

	// helper function to transpose an array
	auto transpose_into = [](double* out_x_t, clapack_int ldx_t, const double* x, clapack_int ldx, clapack_int n, int matrix_layout, char uplo) {
		bool is_colmaj = (LAPACK_COL_MAJOR == matrix_layout);
		bool is_lower = ('L' == uplo);
		clapack_int j_size = minimum(n, ldx_t);

		if(is_colmaj != is_lower) {
			for(clapack_int j = 0; j < j_size; ++j) {
				clapack_int i_size = minimum(j+1, ldx);
				for(clapack_int i = 0; i < i_size; ++i) {
					out_x_t[j+i*ldx_t] = x[i+j*ldx];
				}
			}
		} else {
			clapack_int i_size = minimum(n, ldx);
			for(clapack_int j = 0; j < j_size; ++j) {
				for(clapack_int i = j; i < i_size; ++i) {
					out_x_t[j+i*ldx_t] = x[i+j*ldx];
				}
			}
		}
	};

	// helper function to query optimal working size
	auto query_lwork = [&](clapack_int& out_lwork, clapack_int lda) -> int {
		clapack_int info = 0;
		clapack_int lwork = -1; // flag to query work size
		double work_d; // a length-1 array of working space, as far as `dsytrf_` is concerned
		dsytrf_(&uplo, &n, a, &lda, ipiv, &work_d, &lwork, &info);
		if(info < 0) {
			return info-1;
		} else {
			out_lwork = static_cast<clapack_int>(work_d);
			return 0;
		}
	};

	if(LAPACK_ROW_MAJOR == matrix_layout) {
		clapack_int info = 0;

		// allocate working space
		clapack_int lwork;
		clapack_int lda_t = maximum<clapack_int>(1,n);
		info = query_lwork(lwork, lda_t);
		if(info < 0) return info;
		std::vector<double> work(lwork);

		// call CLAPACK function on the transposed array
		std::vector<double> a_t(lda_t * maximum<clapack_int>(1,n));
		transpose_into(a_t.data(), lda_t, a, lda, n, matrix_layout, uplo);
		dsytrf_(&uplo, &n, a_t.data(), &lda_t, ipiv, work.data(), &lwork, &info);
		if(info < 0) return info-1;
		transpose_into(a, lda, a_t.data(), lda_t, n, LAPACK_COL_MAJOR, uplo);
	} else if(LAPACK_COL_MAJOR == matrix_layout) {
		clapack_int info = 0;

		// allocate working space
		clapack_int lwork;
		info = query_lwork(lwork, lda);
		if(info < 0) return info;
		std::vector<double> work(lwork);

		// call CLAPACK function
		dsytrf_(&uplo, &n, a, &lda, ipiv, work.data(), &lwork, &info);
		if(info < 0) return info-1;
	} else {
		return -1; // argument 1 has an illegal value
	}

	return 0;
}

int LAPACKE_dsytrf(
	int matrix_layout,
	char uplo,
	int n,
	double* a,
	int lda,
	int* ipiv 
) {
	std::vector<clapack_int> cipiv(n);
	int info = LAPACKE_dsytrf(
		matrix_layout,
		uplo,
		static_cast<clapack_int>(n),
		a,
		static_cast<clapack_int>(lda),
		cipiv.data()
	);

	for(int i = 0; i < n; ++i) {
		ipiv[i] = static_cast<int>(cipiv[i]);
	}

	return info;
}

} // namespace clapacke
} // namespace nasoq
