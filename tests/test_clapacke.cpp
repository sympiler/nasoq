#include <catch2/catch.hpp>

#include <Eigen/Dense>
#include <mkl_lapacke.h>
#include <nasoq/clapacke/clapacke.h>


using MatrixXd = Eigen::MatrixXd;
using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


template<typename Derived>
void fill_matrix(Derived& out_matrix) {
	using Int = typename Eigen::DenseIndex;
	for(Int i = 0; i < out_matrix.rows(); ++i) {
		for(Int j = 0; j < out_matrix.cols(); ++j) {
			out_matrix(i,j) = i + j*j;
		}
	}
}

template<typename Derived>
void fill_spd(Derived& out_matrix) {
	using Int = typename Eigen::DenseIndex;
	Int n = out_matrix.rows();
	REQUIRE(n == out_matrix.cols());

	for(Int i = 0; i < n; ++i) {
		for(Int j = 0; j < n; ++j) {
			// out_matrix(i,j) = 1.0 / (i + j + 1); // Hilbert matrix
			out_matrix(i,j) = 1.0 / (2*n - i - j - 1); // Hilbert matrix with reversed indices
		}
	}
}

//----------------------------------------------------------------------------//

TEST_CASE("nasoq::clapacke::LAPACKE_dsytrf") {
	const int n = 4;

	// reference matrix
	MatrixXd A(n, n);
	fill_spd(A);

	// Note that the factorizations aren't required to be unique since the
	// pivoting is allowed to differ, so we compare against MLK my using
	// `LAPACKE_dsytri` to compute the inverse of `A` from the factorization.

	SECTION("col-major") {
		SECTION("upper") {
			MatrixXd A0(n, n);
			fill_spd(A0);
			MatrixXd A1 = A0;
			Eigen::ArrayXi ipiv0(n), ipiv1(n);

			int ret0 = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'U', n, A0.data(), A0.outerStride(), ipiv0.data());
			int ret1 = nasoq::clapacke::LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'U', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK(ret0 == ret1);

			LAPACKE_dsytri(LAPACK_COL_MAJOR, 'U', n, A0.data(), A0.outerStride(), ipiv0.data());
			LAPACKE_dsytri(LAPACK_COL_MAJOR, 'U', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK( A0.isApprox(A1) );
		}

		SECTION("lower") {
			MatrixXd A0(n, n);
			fill_spd(A0);
			MatrixXd A1 = A0;
			Eigen::ArrayXi ipiv0(n), ipiv1(n);

			int ret0 = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, A0.data(), A0.outerStride(), ipiv0.data());
			int ret1 = nasoq::clapacke::LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK(ret0 == ret1);

			LAPACKE_dsytri(LAPACK_COL_MAJOR, 'L', n, A0.data(), A0.outerStride(), ipiv0.data());
			LAPACKE_dsytri(LAPACK_COL_MAJOR, 'L', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK( A0.isApprox(A1) );
		}
	}

	SECTION("row-major") {
		SECTION("upper") {
			RowMatrixXd A0(n, n);
			fill_spd(A0);
			RowMatrixXd A1 = A0;
			Eigen::ArrayXi ipiv0(n), ipiv1(n);

			int ret0 = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, A0.data(), A0.outerStride(), ipiv0.data());
			int ret1 = nasoq::clapacke::LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK(ret0 == ret1);

			LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, A0.data(), A0.outerStride(), ipiv0.data());
			LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK( A0.isApprox(A1) );
		}

		SECTION("lower") {
			RowMatrixXd A0(n, n);
			fill_spd(A0);
			RowMatrixXd A1 = A0;
			Eigen::ArrayXi ipiv0(n), ipiv1(n);

			int ret0 = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'L', n, A0.data(), A0.outerStride(), ipiv0.data());
			int ret1 = nasoq::clapacke::LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'L', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK(ret0 == ret1);

			LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'L', n, A0.data(), A0.outerStride(), ipiv0.data());
			LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'L', n, A1.data(), A1.outerStride(), ipiv1.data());
			CHECK( A0.isApprox(A1) );
		}
	}
}


TEST_CASE("nasoq::clapacke::LAPACKE_dlapmt") {
	const int m = 5;
	const int n = 7;

	Eigen::VectorXi K(n);
	for(int i = 0; i < n; ++i) {
		K[i] = n - i; // 1-based reversed indices
	}

	SECTION("col-major") {
		MatrixXd X0(m, n);
		fill_matrix(X0);
		MatrixXd X1 = X0;

		int ldx = static_cast<int>(X0.outerStride());
		double* x0 = X0.data();
		double* x1 = X1.data();
		int* k = K.data();

		SECTION("forward") {
			int ret0 = LAPACKE_dlapmt(LAPACK_COL_MAJOR, 1, m, n, x0, ldx, k);
			int ret1 = nasoq::clapacke::LAPACKE_dlapmt(LAPACK_COL_MAJOR, 1, m, n, x1, ldx, k);
			CHECK(ret0 == ret1);
			CHECK(X0 == X1);
		}
		SECTION("backward") {
			int ret0 = LAPACKE_dlapmt(LAPACK_COL_MAJOR, 0, m, n, x0, ldx, k);
			int ret1 = nasoq::clapacke::LAPACKE_dlapmt(LAPACK_COL_MAJOR, 0, m, n, x1, ldx, k);
			CHECK(ret0 == ret1);
			CHECK(X0 == X1);
		}
	}

	SECTION("row-major") {
		RowMatrixXd X0(m, n);
		fill_matrix(X0);
		RowMatrixXd X1 = X0;

		int ldx = static_cast<int>(X0.outerStride());
		double* x0 = X0.data();
		double* x1 = X1.data();
		int* k = K.data();

		SECTION("forward") {
			int ret0 = LAPACKE_dlapmt(LAPACK_ROW_MAJOR, 1, m, n, x0, ldx, k);
			int ret1 = nasoq::clapacke::LAPACKE_dlapmt(LAPACK_ROW_MAJOR, 1, m, n, x1, ldx, k);
			CHECK(ret0 == ret1);
			CHECK(X0 == X1);
		}
		SECTION("backward") {
			int ret0 = LAPACKE_dlapmt(LAPACK_ROW_MAJOR, 0, m, n, x0, ldx, k);
			int ret1 = nasoq::clapacke::LAPACKE_dlapmt(LAPACK_ROW_MAJOR, 0, m, n, x1, ldx, k);
			CHECK(ret0 == ret1);
			CHECK(X0 == X1);
		}
	}
}
