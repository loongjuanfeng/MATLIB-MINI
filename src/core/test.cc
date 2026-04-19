#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "matlib.hh"

namespace {

bool nearly_equal(long double a, long double b, long double tol = 1e-6L) {
  return std::fabsl(a - b) <= tol;
}

void assert_true(bool cond, const std::string& message) {
  if (!cond) {
	throw std::runtime_error(message);
  }
}

void assert_near(long double a, long double b, long double tol,
				 const std::string& message) {
  if (!nearly_equal(a, b, tol)) {
	std::ostringstream oss;
	oss << message << " (actual=" << static_cast<double>(a)
		<< ", expected=" << static_cast<double>(b)
		<< ", tol=" << static_cast<double>(tol) << ")";
	throw std::runtime_error(oss.str());
  }
}

void assert_shape(const Matrix& mat, std::size_t r, std::size_t c,
				  const std::string& message) {
  auto sh = mat.shape();
  if (sh.first != r || sh.second != c) {
	std::ostringstream oss;
	oss << message << " (actual=" << sh.first << "x" << sh.second
		<< ", expected=" << r << "x" << c << ")";
	throw std::runtime_error(oss.str());
  }
}

void assert_matrix_near(const Matrix& mat,
						const std::vector<std::vector<long double>>& expected,
						long double tol, const std::string& message) {
  assert_shape(mat, expected.size(), expected.empty() ? 0 : expected[0].size(),
			   message + " shape mismatch");
  for (std::size_t i = 0; i < expected.size(); ++i) {
	for (std::size_t j = 0; j < expected[i].size(); ++j) {
	  assert_near(mat(i, j), expected[i][j], tol, message + " value mismatch");
	}
  }
}

class TestRunner {
 private:
  int passed_ = 0;
  int failed_ = 0;

 public:
  template <typename Fn>
  void run(const std::string& name, Fn&& fn) {
	try {
	  fn();
	  ++passed_;
	  std::cout << "[PASS] " << name << "\n";
	} catch (const std::exception& e) {
	  ++failed_;
	  std::cout << "[FAIL] " << name << " -> " << e.what() << "\n";
	} catch (...) {
	  ++failed_;
	  std::cout << "[FAIL] " << name << " -> unknown error\n";
	}
  }

  int passed() const { return passed_; }
  int failed() const { return failed_; }
};

}  // namespace

int main() {
  TestRunner tr;

  tr.run("Rational arithmetic and comparisons", [] {
	Rational a(1, 2);
	Rational b(1, 3);
	Rational c = a + b;
	assert_true(c == Rational(5, 6), "Rational addition failed");
	assert_true((a - b) == Rational(1, 6), "Rational subtraction failed");
	assert_true((a * b) == Rational(1, 6), "Rational multiplication failed");
	assert_true((a / b) == Rational(3, 2), "Rational division failed");
	assert_true((-a) == Rational(-1, 2), "Rational unary minus failed");

	Rational d = Rational::fromLongDouble(0.25L);
	assert_true(d == Rational(1, 4), "fromLongDouble failed");
	assert_true(!d.isZero(), "isZero failed");
	assert_near(d.toLongDouble(), 0.25L, 1e-12L, "toLongDouble failed");
	assert_true(Rational(1, 2) < Rational(2, 3), "Rational comparison failed");

	Rational x(1, 2);
	x += Rational(1, 2);
	assert_true(x == Rational(1, 1), "operator+= failed");
	x -= Rational(1, 4);
	assert_true(x == Rational(3, 4), "operator-= failed");
	x *= Rational(8, 3);
	assert_true(x == Rational(2, 1), "operator*= failed");
	x /= Rational(2, 1);
	assert_true(x == Rational(1, 1), "operator/= failed");
  });

  tr.run("Constructors and basic metadata", [] {
	Matrix m0;
	assert_shape(m0, 0, 0, "default constructor shape");
	assert_true(m0.dtype() == DType::Fraction, "default dtype");

	Matrix m1(2, 3, 2.0L, 1e-9L, DType::Float);
	assert_true(m1.getRows() == 2 && m1.getCols() == 3,
				"getRows/getCols failed");
	assert_true(m1.rows() == 2 && m1.cols() == 3, "rows/cols failed");
	assert_true(m1.size() == 6, "size failed");
	assert_true(m1.dtype() == DType::Float, "dtype failed");
	assert_true(m1.dtype_str() == "float", "dtype_str failed");
	assert_near(m1.tolerance(), 1e-9L, 1e-15L, "tolerance failed");
	assert_true(Matrix::method_comparison().find("Gaussian elimination") !=
					std::string::npos,
				"method_comparison content");
  });

  tr.run("Element access and casting", [] {
	Matrix A(2, 2, 0.0L, 1e-10L, DType::Fraction);
	A(0, 0) = 1.5L;
	A(0, 1) = Rational(1, 3);
	A.set_value(1, 0, Rational(5, 2));
	A.set_value(1, 1, -4.0L);

	assert_true(A.at_fraction(0, 0) == Rational(3, 2), "at_fraction value");
	assert_near(A.at_float(0, 1), 1.0L / 3.0L, 1e-12L, "at_float value");

	const Matrix& C = A;
	assert_near(C(1, 0), 2.5L, 1e-12L, "const operator() failed");

	auto fd = A.to_float_data();
	assert_near(fd[1][1], -4.0L, 1e-12L, "to_float_data failed");
  });

  tr.run("String output, print and stream operator", [] {
	Matrix A({{1.0L, 2.0L}, {3.0L, 4.0L}}, 1e-10L, DType::Float);
	std::string s = A.str(2);
	assert_true(s.find("1.00") != std::string::npos, "str precision failed");

	std::ostringstream oss;
	oss << A;
	assert_true(oss.str().find("[") != std::string::npos, "operator<< failed");

	A.print(2);
  });

  tr.run("copy reshape transpose T identity I", [] {
	Matrix A({{1.0L, 2.0L}, {3.0L, 4.0L}}, 1e-10L, DType::Float);
	Matrix B = A.copy();
	B(0, 0) = 10.0L;
	assert_near(A(0, 0), 1.0L, 1e-12L, "copy should be deep copy");

	Matrix C = A.reshape({1, 4});
	assert_matrix_near(C, {{1.0L, 2.0L, 3.0L, 4.0L}}, 1e-12L, "reshape(pair)");
	Matrix D = A.reshape(4, 1);
	assert_matrix_near(D, {{1.0L}, {2.0L}, {3.0L}, {4.0L}}, 1e-12L,
					   "reshape(rows,cols)");

	Matrix At = A.transpose();
	Matrix AT = A.T();
	assert_matrix_near(At, {{1.0L, 3.0L}, {2.0L, 4.0L}}, 1e-12L, "transpose");
	assert_matrix_near(AT, {{1.0L, 3.0L}, {2.0L, 4.0L}}, 1e-12L, "T alias");

	Matrix I1 = Matrix::identity(3, 1e-10L, DType::Float);
	Matrix I2 = I(3, 1e-10L, "float");
	assert_shape(I1, 3, 3, "identity shape");
	assert_shape(I2, 3, 3, "I shape");
	assert_near(I1(2, 2), 1.0L, 1e-12L, "identity diagonal");
	assert_near(I2(1, 1), 1.0L, 1e-12L, "I diagonal");
  });

  tr.run("Symmetry and sums", [] {
	Matrix S({{1.0L, 2.0L, 3.0L}, {4.0L, 5.0L, 6.0L}}, 1e-10L, DType::Float);
	assert_true(!S._is_symmetric(), "_is_symmetric false case");
	Matrix Sym({{2.0L, 1.0L}, {1.0L, 2.0L}}, 1e-10L, DType::Float);
	assert_true(Sym.is_symmetric(), "is_symmetric true case");

	assert_matrix_near(S.sum(), {{21.0L}}, 1e-12L, "sum none");
	assert_matrix_near(S.sum(0), {{5.0L, 7.0L, 9.0L}}, 1e-12L, "sum axis0");
	assert_matrix_near(S.sum(1), {{6.0L}, {15.0L}}, 1e-12L, "sum axis1");
  });

  tr.run("Concatenate and Kronecker product", [] {
	Matrix A({{1.0L, 2.0L}, {3.0L, 4.0L}}, 1e-10L, DType::Float);
	Matrix B({{5.0L, 6.0L}}, 1e-10L, DType::Float);
	Matrix C = A.concatenate(B, 0);
	assert_shape(C, 3, 2, "row concatenate shape");

	Matrix D({{7.0L}, {8.0L}}, 1e-10L, DType::Float);
	Matrix E = A.concatenate(D, 1);
	assert_shape(E, 2, 3, "col concatenate shape");

	Matrix G = concatenate({A, B}, 0);
	assert_shape(G, 3, 2, "global concatenate shape");

	Matrix K1({{1.0L, 2.0L}, {3.0L, 4.0L}}, 1e-10L, DType::Float);
	Matrix K2({{0.0L, 5.0L}, {6.0L, 7.0L}}, 1e-10L, DType::Float);
	Matrix KP = K1.kronecker_product(K2);
	Matrix KP2 = K1.Kronecker_product(K2);
	assert_shape(KP, 4, 4, "kronecker shape");
	assert_near(KP(0, 1), 5.0L, 1e-12L, "kronecker value");
	assert_near(KP(3, 3), 28.0L, 1e-12L, "kronecker value");
	assert_near(KP2(3, 3), 28.0L, 1e-12L, "Kronecker alias");
  });

  tr.run("Scalar and matrix arithmetic", [] {
	Matrix A({{1.0L, 2.0L}, {3.0L, 4.0L}}, 1e-10L, DType::Float);
	Matrix B = A.num_product(2LL);
	Matrix C = A.num_product(2.5L);
	Matrix D = A * 2LL;
	Matrix E = A * 2.5L;
	Matrix F = 3LL * A;
	Matrix G = 1.5L * A;
	assert_near(B(1, 1), 8.0L, 1e-12L, "num_product int");
	assert_near(C(1, 1), 10.0L, 1e-12L, "num_product float");
	assert_near(D(0, 1), 4.0L, 1e-12L, "operator* int");
	assert_near(E(0, 1), 5.0L, 1e-12L, "operator* long double");
	assert_near(F(1, 0), 9.0L, 1e-12L, "friend operator* int");
	assert_near(G(1, 0), 4.5L, 1e-12L, "friend operator* long double");

	Matrix H({{1.0L, 1.0L}, {1.0L, 1.0L}}, 1e-10L, DType::Float);
	Matrix add = A + H;
	Matrix sub = A - H;
	assert_near(add(1, 1), 5.0L, 1e-12L, "operator+ failed");
	assert_near(sub(1, 1), 3.0L, 1e-12L, "operator- failed");

	Matrix mul = A * H;
	assert_matrix_near(mul, {{3.0L, 3.0L}, {7.0L, 7.0L}}, 1e-12L,
					   "operator*(matrix)");
  });

  tr.run("Power functions", [] {
	Matrix A({{2.0L, 0.0L}, {0.0L, 3.0L}}, 1e-10L, DType::Float);
	Matrix p3 = A.power(3);
	Matrix p2 = A.pow(2);
	assert_matrix_near(p3, {{8.0L, 0.0L}, {0.0L, 27.0L}}, 1e-12L, "power");
	assert_matrix_near(p2, {{4.0L, 0.0L}, {0.0L, 9.0L}}, 1e-12L, "pow alias");
  });

  tr.run("Gaussian elimination determinant det rank inverse", [] {
	Matrix A({{2.0L, 1.0L}, {4.0L, 3.0L}}, 1e-10L, DType::Float);
	auto ge = A.gauss_elimination(true);
	assert_shape(ge.first, 2, 2, "gauss_elimination shape");
	assert_true(ge.second >= 0, "swap_count non-negative");

	long double detA = A.determinant();
	assert_near(detA, 2.0L, 1e-9L, "determinant failed");
	assert_near(A.det(), 2.0L, 1e-9L, "det alias failed");
	assert_true(A.rank() == 2, "rank full failed");

	Matrix invA = A.inverse();
	Matrix I2 = A * invA;
	assert_matrix_near(I2, {{1.0L, 0.0L}, {0.0L, 1.0L}}, 1e-6L,
					   "inverse reconstruction");

	Matrix B({{1.0L, 2.0L}, {2.0L, 4.0L}}, 1e-10L, DType::Float);
	assert_true(B.rank() == 1, "rank singular failed");
  });

  tr.run("LU decompositions and wrappers", [] {
	Matrix A({{4.0L, 3.0L}, {6.0L, 3.0L}}, 1e-10L, DType::Float);
	auto plu = A.lu_decomposition(true);
	Matrix P = std::get<0>(plu);
	Matrix L = std::get<1>(plu);
	Matrix U = std::get<2>(plu);
	assert_matrix_near(P * A, (L * U).to_float_data(), 1e-6L, "PLU relation");

	auto lu = A.lu_decomposition_no_pivot();
	Matrix L2 = lu.first;
	Matrix U2 = lu.second;
	assert_matrix_near(L2 * U2, A.to_float_data(), 1e-6L,
					   "LU no pivot relation");

	Matrix L3, U3;
	A.luDecomposition(L3, U3);
	assert_matrix_near(L3 * U3, A.to_float_data(), 1e-6L,
					   "luDecomposition wrapper");
  });

  tr.run("Cholesky decomposition", [] {
	Matrix A({{4.0L, 2.0L}, {2.0L, 3.0L}}, 1e-10L, DType::Float);
	auto ch = A.cholesky_decomposition();
	Matrix L = ch.first;
	Matrix LT = ch.second;
	assert_matrix_near(L * LT, A.to_float_data(), 1e-6L,
					   "cholesky reconstruction");
  });

  tr.run("Gram-Schmidt and QR families", [] {
	Matrix A({{1.0L, 1.0L}, {1.0L, -1.0L}, {1.0L, 0.0L}}, 1e-10L, DType::Float);
	Matrix Q1 = A.gram_schmidt("mgs");
	Matrix Q2 = A.gram_schmidt("cgs");
	Matrix I1 = Q1.transpose() * Q1;
	Matrix I2 = Q2.transpose() * Q2;
	assert_matrix_near(I1, {{1.0L, 0.0L}, {0.0L, 1.0L}}, 1e-5L,
					   "MGS orthonormality");
	assert_matrix_near(I2, {{1.0L, 0.0L}, {0.0L, 1.0L}}, 1e-5L,
					   "CGS orthonormality");

	auto hh = A._qr_householder();
	assert_matrix_near(hh.first * hh.second, A.to_float_data(), 1e-5L,
					   "_qr_householder reconstruction");

	auto qr = A.qr_decomposition("householder");
	assert_matrix_near(qr.first * qr.second, A.to_float_data(), 1e-5L,
					   "qr_decomposition reconstruction");

	auto qr2 = A.QR_decomposition("mgs");
	assert_matrix_near(qr2.first * qr2.second, A.to_float_data(), 1e-5L,
					   "QR_decomposition alias reconstruction");

	Matrix Q, R;
	A.qrDecomposition(Q, R);
	assert_matrix_near(Q * R, A.to_float_data(), 1e-5L,
					   "qrDecomposition wrapper reconstruction");
  });

  tr.run("Eigenvalues and eigenpairs", [] {
	Matrix A({{2.0L, 1.0L}, {1.0L, 2.0L}}, 1e-10L, DType::Float);

	auto vals_qr = A.eigenvalues_qr(120);
	assert_true(vals_qr.size() == 2, "eigenvalues_qr size");
	std::sort(vals_qr.begin(), vals_qr.end(), std::greater<long double>());
	assert_near(vals_qr[0], 3.0L, 1e-3L, "eigenvalues_qr lambda1");
	assert_near(vals_qr[1], 1.0L, 1e-3L, "eigenvalues_qr lambda2");

	auto jac = A.eigenpairs_jacobi(2000, true);
	auto vals_j = jac.first;
	Matrix V = jac.second;
	assert_near(vals_j[0], 3.0L, 1e-6L, "jacobi lambda1");
	assert_near(vals_j[1], 1.0L, 1e-6L, "jacobi lambda2");

	std::vector<long double> v0 = {V(0, 0), V(1, 0)};
	Matrix vec({{v0[0]}, {v0[1]}}, 1e-10L, DType::Float);
	Matrix Av = A * vec;
	Matrix lv = vec * vals_j[0];
	assert_matrix_near(Av, lv.to_float_data(), 1e-5L,
					   "jacobi eigenvector relation");

	auto auto_pair = A.eigenpairs("auto");
	assert_true(auto_pair.first.size() == 2, "eigenpairs auto values");
	assert_shape(auto_pair.second, 2, 2, "eigenpairs auto vectors shape");

	auto qr_pair = A.eigenpairs("qr", 60);
	assert_true(qr_pair.first.size() == 2, "eigenpairs qr values");
	assert_shape(qr_pair.second, 0, 0, "eigenpairs qr vectors empty");

	auto vals = A.eigenvalues("auto");
	assert_true(vals.size() == 2, "eigenvalues size");
  });

  tr.run("SVD and alias svd", [] {
	Matrix A({{3.0L, 1.0L}, {1.0L, 3.0L}}, 1e-10L, DType::Float);
	auto s1 = A.singular_value_decomposition();
	Matrix U = std::get<0>(s1);
	Matrix S = std::get<1>(s1);
	Matrix VT = std::get<2>(s1);
	assert_shape(U, 2, 2, "SVD U shape");
	assert_shape(S, 2, 2, "SVD Sigma shape");
	assert_shape(VT, 2, 2, "SVD VT shape");
	assert_matrix_near(U * S * VT, A.to_float_data(), 1e-5L,
					   "SVD reconstruction");

	auto s2 = A.svd();
	Matrix U2 = std::get<0>(s2);
	Matrix S2 = std::get<1>(s2);
	Matrix VT2 = std::get<2>(s2);
	assert_matrix_near(U2 * S2 * VT2, A.to_float_data(), 1e-5L,
					   "svd alias reconstruction");
  });

  tr.run("Jordan normal form and alias", [] {
	Matrix A({{2.0L, 1.0L}, {0.0L, 2.0L}}, 1e-10L, DType::Fraction);
	Matrix J1 = A.jordan_normal_form();
	Matrix J2 = A.jordanNormalForm();
	assert_matrix_near(J1, {{2.0L, 1.0L}, {0.0L, 2.0L}}, 1e-9L,
					   "jordan_normal_form");
	assert_matrix_near(J2, {{2.0L, 1.0L}, {0.0L, 2.0L}}, 1e-9L,
					   "jordanNormalForm alias");
  });

  tr.run("dominantEigen and fullRankDecomposition", [] {
	Matrix A({{2.0L, 0.0L}, {0.0L, 1.0L}}, 1e-10L, DType::Float);
	long double lambda = 0.0L;
	Matrix vec;
	A.dominantEigen(lambda, vec, 80);
	assert_near(lambda, 2.0L, 1e-4L, "dominantEigen value");

	long double nrm = std::sqrt(vec(0, 0) * vec(0, 0) + vec(1, 0) * vec(1, 0));
	assert_near(nrm, 1.0L, 1e-4L, "dominantEigen normalized vector");

	Matrix B({{1.0L, 2.0L, 3.0L}, {2.0L, 4.0L, 6.0L}, {1.0L, 1.0L, 1.0L}},
			 1e-10L, DType::Float);
	Matrix F, G;
	B.fullRankDecomposition(F, G);
	assert_matrix_near(F * G, B.to_float_data(), 1e-5L,
					   "fullRankDecomposition reconstruction");
  });

  std::cout << "\n========== Test Summary ==========" << "\n";
  std::cout << "Passed: " << tr.passed() << "\n";
  std::cout << "Failed: " << tr.failed() << "\n";
  std::cout << "==================================" << "\n";

  return tr.failed() == 0 ? 0 : 1;
}
