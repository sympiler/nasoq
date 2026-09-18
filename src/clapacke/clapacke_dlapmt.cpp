#include "nasoq/clapacke/clapacke.h"

#include <vector>

#if defined(ACCELERATE)
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include "f2c.h"
#include "clapack.h"
}
#endif


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
	clapack_logical forwrd,
	clapack_int m,
	clapack_int n,
	double* x,
	clapack_int ldx,
	clapack_int* k
) {
	auto transpose_into = [](double* out_x_t, clapack_int ldx_t, const double* x, clapack_int ldx, clapack_int m, clapack_int n, int matrix_layout) {
		clapack_int i_max, j_max;
		if(LAPACK_COL_MAJOR == matrix_layout) {
			j_max = n;
			i_max = m;
		} else {
			j_max = m;
			i_max = n;
		}
		j_max = minimum(j_max, ldx_t);
		i_max = minimum(i_max, ldx);

		for(clapack_int i = 0; i < i_max; ++i) {
			for(clapack_int j = 0; j < j_max; ++j) {
				out_x_t[i*ldx_t + j] = x[j*ldx + i];
			}
		}
	};

	if(LAPACK_ROW_MAJOR == matrix_layout) {
		clapack_int ldx_t = maximum<clapack_int>(1,m);
		if(ldx < n) return -6;
		std::vector<double> x_t(ldx_t * maximum<clapack_int>(1,n));
		transpose_into(x_t.data(), ldx_t, x, ldx, m, n, matrix_layout);
#if defined(ACCELERATE)
		// Accelerate's `dlapmt_` matches reference LAPACK exactly: it is a
		// `void`-returning subroutine with no INFO output argument (unlike
		// CLAPACK's f2c-translated version, which synthesizes an `integer`
		// return value), so there is nothing to check here.
		dlapmt_(&forwrd, &m, &n, x_t.data(), &ldx_t, k);
#else
		clapack_int info = dlapmt_(&forwrd, &m, &n, x_t.data(), &ldx_t, k);
		if(info < 0) return info-1;
#endif
		transpose_into(x, ldx, x_t.data(), ldx_t, m, n, LAPACK_COL_MAJOR);
	} else if(LAPACK_COL_MAJOR == matrix_layout) {
#if defined(ACCELERATE)
		dlapmt_(&forwrd, &m, &n, x, &ldx, k);
#else
		clapack_int info = dlapmt_(&forwrd, &m, &n, x, &ldx, k);
		if(info < 0) return info-1;
#endif
	} else {
		return -1; // argument 1 has an illegal value
	}
	return 0;
}

// Note: the `int`-typed overload of LAPACKE_dlapmt declared in clapacke.h
// is intentionally not redefined here: now that `clapack_int`/
// `clapack_logical` are `int` (see clapacke.h), it would be an
// exact-signature redefinition of the clapack_int overload above.

} // namespace clapacke
} // namespace nasoq
