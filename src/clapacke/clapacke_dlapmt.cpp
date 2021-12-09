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
		clapack_int info = dlapmt_(&forwrd, &m, &n, x_t.data(), &ldx_t, k);
		if(info < 0) return info-1;
		transpose_into(x, ldx, x_t.data(), ldx_t, m, n, LAPACK_COL_MAJOR);
	} else if(LAPACK_COL_MAJOR == matrix_layout) {
		clapack_int info = dlapmt_(&forwrd, &m, &n, x, &ldx, k);
		if(info < 0) return info-1;
	} else {
		return -1; // argument 1 has an illegal value
	}
	return 0;
}

int LAPACKE_dlapmt(
	int matrix_layout,
	int forwrd,
	int m,
	int n,
	double* x,
	int ldx,
	int* k
) {
	std::vector<clapack_int> ck(n);
	for(int i = 0; i < n; ++i) {
		ck[i] = static_cast<clapack_int>(k[i]);
	}
	
	int info = LAPACKE_dlapmt(
		matrix_layout,
		forwrd,
		static_cast<clapack_int>(m),
		static_cast<clapack_int>(n),
		x,
		static_cast<clapack_int>(ldx),
		ck.data()
	);

	return info;
}

} // namespace clapacke
} // namespace nasoq
