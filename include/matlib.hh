#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class Rational {
private:
	std::int64_t num_;
	std::int64_t den_;

	static std::int64_t gcd64(std::int64_t a, std::int64_t b) {
		a = std::llabs(a);
		b = std::llabs(b);
		while (b != 0) {
			std::int64_t t = a % b;
			a = b;
			b = t;
		}
		return a == 0 ? 1 : a;
	}

	void normalize() {
		if (den_ == 0) {
			throw std::invalid_argument("Rational denominator cannot be zero.");
		}
		if (den_ < 0) {
			num_ = -num_;
			den_ = -den_;
		}
		std::int64_t g = gcd64(num_, den_);
		num_ /= g;
		den_ /= g;
	}

public:
	Rational(std::int64_t numerator = 0, std::int64_t denominator = 1)
		: num_(numerator), den_(denominator) {
		normalize();
	}

	static Rational fromLongDouble(long double value, std::int64_t max_den = 1000000) {
		if (!std::isfinite(static_cast<double>(value))) {
			throw std::invalid_argument("Cannot cast non-finite value to Rational.");
		}
		long double rounded = std::round(value);
		if (std::fabs(value - rounded) <= 1e-12L) {
			return Rational(static_cast<std::int64_t>(rounded), 1);
		}

		long double scaled = std::round(value * static_cast<long double>(max_den));
		if (scaled > static_cast<long double>(std::numeric_limits<std::int64_t>::max()) ||
			scaled < static_cast<long double>(std::numeric_limits<std::int64_t>::min())) {
			throw std::overflow_error("Rational cast overflow.");
		}
		return Rational(static_cast<std::int64_t>(scaled), max_den);
	}

	std::int64_t numerator() const { return num_; }
	std::int64_t denominator() const { return den_; }

	bool isZero() const { return num_ == 0; }

	long double toLongDouble() const {
		return static_cast<long double>(num_) / static_cast<long double>(den_);
	}

	Rational operator-() const {
		return Rational(-num_, den_);
	}

	friend Rational operator+(const Rational& a, const Rational& b) {
		std::int64_t num = a.num_ * b.den_ + b.num_ * a.den_;
		std::int64_t den = a.den_ * b.den_;
		return Rational(num, den);
	}

	friend Rational operator-(const Rational& a, const Rational& b) {
		std::int64_t num = a.num_ * b.den_ - b.num_ * a.den_;
		std::int64_t den = a.den_ * b.den_;
		return Rational(num, den);
	}

	friend Rational operator*(const Rational& a, const Rational& b) {
		std::int64_t num = a.num_ * b.num_;
		std::int64_t den = a.den_ * b.den_;
		return Rational(num, den);
	}

	friend Rational operator/(const Rational& a, const Rational& b) {
		if (b.num_ == 0) {
			throw std::invalid_argument("Division by zero Rational.");
		}
		std::int64_t num = a.num_ * b.den_;
		std::int64_t den = a.den_ * b.num_;
		return Rational(num, den);
	}

	Rational& operator+=(const Rational& other) {
		*this = *this + other;
		return *this;
	}

	Rational& operator-=(const Rational& other) {
		*this = *this - other;
		return *this;
	}

	Rational& operator*=(const Rational& other) {
		*this = *this * other;
		return *this;
	}

	Rational& operator/=(const Rational& other) {
		*this = *this / other;
		return *this;
	}

	bool operator==(const Rational& other) const {
		return num_ == other.num_ && den_ == other.den_;
	}

	bool operator!=(const Rational& other) const {
		return !(*this == other);
	}

	bool operator<(const Rational& other) const {
		return toLongDouble() < other.toLongDouble();
	}
};

enum class DType {
	Fraction,
	Float,
};

static std::string to_lower(std::string text) {
	for (char& ch : text) {
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	return text;
}

static DType parse_dtype(const std::string& dtype) {
	std::string key = to_lower(dtype);
	if (key == "fraction") {
		return DType::Fraction;
	}
	if (key == "float") {
		return DType::Float;
	}
	throw std::invalid_argument("dtype must be 'fraction' or 'float'.");
}

static std::string dtype_to_string(DType dtype) {
	return dtype == DType::Fraction ? "fraction" : "float";
}

class Matrix {
public:
	class CellProxy {
	private:
		Matrix& owner_;
		std::size_t row_;
		std::size_t col_;

	public:
		CellProxy(Matrix& owner, std::size_t row, std::size_t col)
			: owner_(owner), row_(row), col_(col) {}

		CellProxy& operator=(long double value) {
			owner_.set_value(row_, col_, value);
			return *this;
		}

		CellProxy& operator=(const Rational& value) {
			owner_.set_value(row_, col_, value);
			return *this;
		}

		operator long double() const {
			return owner_.at_float(row_, col_);
		}

		operator Rational() const {
			return owner_.at_fraction(row_, col_);
		}
	};

private:
	std::size_t rows_;
	std::size_t cols_;
	DType dtype_;
	long double tolerance_;
	std::vector<std::vector<Rational>> frac_data_;
	std::vector<std::vector<long double>> float_data_;

	static bool near_zero(long double value, long double tol) {
		return std::fabs(value) <= tol;
	}

	static bool near_zero(const Rational& value, long double /*tol*/) {
		return value.isZero();
	}

	static long double dot(const std::vector<long double>& v1, const std::vector<long double>& v2) {
		if (v1.size() != v2.size()) {
			throw std::invalid_argument("Dot product vectors must have the same length.");
		}
		long double out = 0.0L;
		for (std::size_t i = 0; i < v1.size(); ++i) {
			out += v1[i] * v2[i];
		}
		return out;
	}

	static long double norm(const std::vector<long double>& vec) {
		long double n = dot(vec, vec);
		return std::sqrt(std::max(0.0L, n));
	}

	static std::vector<long double> normalize_vec(const std::vector<long double>& vec, long double tol) {
		long double nrm = norm(vec);
		if (nrm <= tol) {
			throw std::invalid_argument("Cannot normalize a near-zero vector.");
		}
		std::vector<long double> out(vec.size(), 0.0L);
		for (std::size_t i = 0; i < vec.size(); ++i) {
			out[i] = vec[i] / nrm;
		}
		return out;
	}

	static std::vector<long double> mat_vec_mul(
		const std::vector<std::vector<long double>>& mat,
		const std::vector<long double>& vec
	) {
		if (mat.empty()) {
			return {};
		}
		if (mat[0].size() != vec.size()) {
			throw std::invalid_argument("Matrix-vector shape mismatch.");
		}
		std::vector<long double> out(mat.size(), 0.0L);
		for (std::size_t i = 0; i < mat.size(); ++i) {
			for (std::size_t j = 0; j < vec.size(); ++j) {
				out[i] += mat[i][j] * vec[j];
			}
		}
		return out;
	}

	static std::vector<std::vector<long double>> complete_orthonormal_basis(
		const std::vector<std::vector<long double>>& base_vectors,
		std::size_t dimension,
		long double tol
	) {
		std::vector<std::vector<long double>> basis;

		for (const auto& vec : base_vectors) {
			if (vec.size() != dimension) {
				throw std::invalid_argument("Base vectors dimension mismatch.");
			}
			std::vector<long double> v = vec;
			for (const auto& q : basis) {
				long double coeff = dot(v, q);
				for (std::size_t i = 0; i < dimension; ++i) {
					v[i] -= coeff * q[i];
				}
			}
			long double nrm = norm(v);
			if (nrm > tol) {
				for (long double& x : v) {
					x /= nrm;
				}
				basis.push_back(v);
			}
		}

		for (std::size_t idx = 0; idx < dimension; ++idx) {
			std::vector<long double> v(dimension, 0.0L);
			v[idx] = 1.0L;
			for (const auto& q : basis) {
				long double coeff = dot(v, q);
				for (std::size_t i = 0; i < dimension; ++i) {
					v[i] -= coeff * q[i];
				}
			}
			long double nrm = norm(v);
			if (nrm > tol) {
				for (long double& x : v) {
					x /= nrm;
				}
				basis.push_back(v);
			}
			if (basis.size() == dimension) {
				break;
			}
		}

		if (basis.size() < dimension) {
			throw std::runtime_error("Failed to complete an orthonormal basis.");
		}
		return basis;
	}

	static std::vector<std::int64_t> get_factors(std::int64_t num) {
		num = std::llabs(num);
		if (num == 0) {
			return {1};
		}
		std::vector<std::int64_t> out;
		for (std::int64_t i = 1; i * i <= num; ++i) {
			if (num % i == 0) {
				out.push_back(i);
				if (i != num / i) {
					out.push_back(num / i);
				}
			}
		}
		return out;
	}

	static Rational evaluate_poly(const std::vector<Rational>& poly, const Rational& x) {
		Rational res(0);
		for (std::size_t idx = poly.size(); idx-- > 0;) {
			res = res * x + poly[idx];
		}
		return res;
	}

	static std::vector<Rational> synthetic_div(const std::vector<Rational>& poly, const Rational& root) {
		std::vector<Rational> new_poly;
		Rational val(0);
		for (std::size_t idx = poly.size(); idx-- > 0;) {
			val = val * root + poly[idx];
			new_poly.push_back(val);
		}
		if (!new_poly.empty()) {
			new_poly.pop_back();
		}
		std::reverse(new_poly.begin(), new_poly.end());
		return new_poly;
	}

	void validate_indices(std::size_t row, std::size_t col) const {
		if (row >= rows_ || col >= cols_) {
			throw std::out_of_range("Matrix index out of range.");
		}
	}

	void swap_rows(std::size_t r1, std::size_t r2) {
		if (r1 == r2) {
			return;
		}
		if (dtype_ == DType::Fraction) {
			std::swap(frac_data_[r1], frac_data_[r2]);
		} else {
			std::swap(float_data_[r1], float_data_[r2]);
		}
	}

	Matrix to_fraction_matrix() const {
		if (dtype_ == DType::Fraction) {
			return *this;
		}
		std::vector<std::vector<Rational>> data(rows_, std::vector<Rational>(cols_, Rational(0)));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				data[i][j] = Rational::fromLongDouble(float_data_[i][j]);
			}
		}
		return Matrix(data, tolerance_);
	}

public:
	Matrix()
		: rows_(0), cols_(0), dtype_(DType::Fraction), tolerance_(1e-10L) {}

	Matrix(
		std::size_t rows,
		std::size_t cols,
		long double init_value = 0.0L,
		long double tolerance = 1e-10L,
		DType dtype = DType::Fraction
	)
		: rows_(rows), cols_(cols), dtype_(dtype), tolerance_(tolerance) {
		if (dtype_ == DType::Fraction) {
			Rational fill = Rational::fromLongDouble(init_value);
			frac_data_.assign(rows_, std::vector<Rational>(cols_, fill));
		} else {
			float_data_.assign(rows_, std::vector<long double>(cols_, init_value));
		}
	}

	Matrix(
		const std::vector<std::vector<long double>>& data,
		long double tolerance = 1e-10L,
		DType dtype = DType::Fraction
	)
		: rows_(data.size()),
		  cols_(data.empty() ? 0 : data.front().size()),
		  dtype_(dtype),
		  tolerance_(tolerance) {
		for (const auto& row : data) {
			if (row.size() != cols_) {
				throw std::invalid_argument("All rows in data must have the same length.");
			}
		}

		if (dtype_ == DType::Fraction) {
			frac_data_.assign(rows_, std::vector<Rational>(cols_, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					frac_data_[i][j] = Rational::fromLongDouble(data[i][j]);
				}
			}
		} else {
			float_data_ = data;
		}
	}

	Matrix(
		const std::vector<std::vector<Rational>>& data,
		long double tolerance = 1e-10L
	)
		: rows_(data.size()),
		  cols_(data.empty() ? 0 : data.front().size()),
		  dtype_(DType::Fraction),
		  tolerance_(tolerance),
		  frac_data_(data) {
		for (const auto& row : data) {
			if (row.size() != cols_) {
				throw std::invalid_argument("All rows in data must have the same length.");
			}
		}
	}

	int getRows() const { return static_cast<int>(rows_); }
	int getCols() const { return static_cast<int>(cols_); }
	std::size_t rows() const { return rows_; }
	std::size_t cols() const { return cols_; }
	std::pair<std::size_t, std::size_t> shape() const { return {rows_, cols_}; }
	std::size_t size() const { return rows_ * cols_; }
	long double tolerance() const { return tolerance_; }
	DType dtype() const { return dtype_; }
	std::string dtype_str() const { return dtype_to_string(dtype_); }

	long double at_float(std::size_t row, std::size_t col) const {
		validate_indices(row, col);
		if (dtype_ == DType::Fraction) {
			return frac_data_[row][col].toLongDouble();
		}
		return float_data_[row][col];
	}

	Rational at_fraction(std::size_t row, std::size_t col) const {
		validate_indices(row, col);
		if (dtype_ == DType::Fraction) {
			return frac_data_[row][col];
		}
		return Rational::fromLongDouble(float_data_[row][col]);
	}

	void set_value(std::size_t row, std::size_t col, long double value) {
		validate_indices(row, col);
		if (dtype_ == DType::Fraction) {
			frac_data_[row][col] = Rational::fromLongDouble(value);
		} else {
			float_data_[row][col] = value;
		}
	}

	void set_value(std::size_t row, std::size_t col, const Rational& value) {
		validate_indices(row, col);
		if (dtype_ == DType::Fraction) {
			frac_data_[row][col] = value;
		} else {
			float_data_[row][col] = value.toLongDouble();
		}
	}

	CellProxy operator()(std::size_t row, std::size_t col) {
		validate_indices(row, col);
		return CellProxy(*this, row, col);
	}

	long double operator()(std::size_t row, std::size_t col) const {
		return at_float(row, col);
	}

	std::vector<std::vector<long double>> to_float_data() const {
		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(cols_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[i][j] = at_float(i, j);
			}
		}
		return out;
	}

	static std::string method_comparison() {
		return
			"# Method Comparison\n"
			"- Gaussian elimination: partial pivoting is more stable than no-pivot elimination; use pivoting by default.\n"
			"- QR decomposition: Modified Gram-Schmidt is usually better than Classical Gram-Schmidt in finite precision.\n"
			"- QR with Householder is typically the most numerically stable general-purpose QR.\n"
			"- Eigenvalues: Jacobi is preferred for real symmetric matrices (stable, gives orthonormal eigenvectors).\n"
			"- Eigenvalues: QR iteration is more general for non-symmetric matrices, but eigenvectors are not returned here.\n"
			"- SVD: using A^T A for V and then U = A v_i / sigma_i is more robust than independently matching U from A A^T.\n"
			"- Cholesky is best for symmetric positive-definite matrices; LU is more general.";
	}

	bool _is_symmetric() const {
		if (rows_ != cols_) {
			return false;
		}
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = i + 1; j < cols_; ++j) {
				if (std::fabs(at_float(i, j) - at_float(j, i)) > tolerance_) {
					return false;
				}
			}
		}
		return true;
	}

	bool is_symmetric() const {
		return _is_symmetric();
	}

	Matrix copy() const {
		return *this;
	}

	Matrix reshape(const std::pair<std::size_t, std::size_t>& newdim) const {
		std::size_t nr = newdim.first;
		std::size_t nc = newdim.second;
		if (rows_ * cols_ != nr * nc) {
			throw std::invalid_argument("Cannot reshape to a different total number of elements.");
		}

		if (dtype_ == DType::Fraction) {
			std::vector<Rational> flat;
			flat.reserve(rows_ * cols_);
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					flat.push_back(frac_data_[i][j]);
				}
			}

			std::vector<std::vector<Rational>> out(nr, std::vector<Rational>(nc, Rational(0)));
			std::size_t idx = 0;
			for (std::size_t i = 0; i < nr; ++i) {
				for (std::size_t j = 0; j < nc; ++j) {
					out[i][j] = flat[idx++];
				}
			}
			return Matrix(out, tolerance_);
		}

		std::vector<long double> flat;
		flat.reserve(rows_ * cols_);
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				flat.push_back(float_data_[i][j]);
			}
		}

		std::vector<std::vector<long double>> out(nr, std::vector<long double>(nc, 0.0L));
		std::size_t idx = 0;
		for (std::size_t i = 0; i < nr; ++i) {
			for (std::size_t j = 0; j < nc; ++j) {
				out[i][j] = flat[idx++];
			}
		}
		return Matrix(out, tolerance_, DType::Float);
	}

	Matrix reshape(std::size_t rows, std::size_t cols) const {
		return reshape({rows, cols});
	}

	Matrix transpose() const {
		if (dtype_ == DType::Fraction) {
			std::vector<std::vector<Rational>> out(cols_, std::vector<Rational>(rows_, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					out[j][i] = frac_data_[i][j];
				}
			}
			return Matrix(out, tolerance_);
		}

		std::vector<std::vector<long double>> out(cols_, std::vector<long double>(rows_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[j][i] = float_data_[i][j];
			}
		}
		return Matrix(out, tolerance_, DType::Float);
	}

	Matrix T() const {
		return transpose();
	}

	static Matrix identity(
		std::size_t n,
		long double tolerance = 1e-10L,
		DType dtype = DType::Fraction
	) {
		if (dtype == DType::Fraction) {
			std::vector<std::vector<Rational>> out(n, std::vector<Rational>(n, Rational(0)));
			for (std::size_t i = 0; i < n; ++i) {
				out[i][i] = Rational(1);
			}
			return Matrix(out, tolerance);
		}

		std::vector<std::vector<long double>> out(n, std::vector<long double>(n, 0.0L));
		for (std::size_t i = 0; i < n; ++i) {
			out[i][i] = 1.0L;
		}
		return Matrix(out, tolerance, DType::Float);
	}

	Matrix sum(const std::optional<int>& axis = std::nullopt) const {
		if (axis.has_value() && axis.value() != 0 && axis.value() != 1) {
			throw std::invalid_argument("axis must be None, 0, or 1.");
		}

		if (!axis.has_value()) {
			if (dtype_ == DType::Fraction) {
				Rational total(0);
				for (std::size_t i = 0; i < rows_; ++i) {
					for (std::size_t j = 0; j < cols_; ++j) {
						total += frac_data_[i][j];
					}
				}
				return Matrix(std::vector<std::vector<Rational>>{{total}}, tolerance_);
			}

			long double total = 0.0L;
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					total += float_data_[i][j];
				}
			}
			return Matrix(std::vector<std::vector<long double>>{{total}}, tolerance_, DType::Float);
		}

		if (axis.value() == 0) {
			if (dtype_ == DType::Fraction) {
				std::vector<std::vector<Rational>> out(1, std::vector<Rational>(cols_, Rational(0)));
				for (std::size_t j = 0; j < cols_; ++j) {
					Rational col_sum(0);
					for (std::size_t i = 0; i < rows_; ++i) {
						col_sum += frac_data_[i][j];
					}
					out[0][j] = col_sum;
				}
				return Matrix(out, tolerance_);
			}

			std::vector<std::vector<long double>> out(1, std::vector<long double>(cols_, 0.0L));
			for (std::size_t j = 0; j < cols_; ++j) {
				for (std::size_t i = 0; i < rows_; ++i) {
					out[0][j] += float_data_[i][j];
				}
			}
			return Matrix(out, tolerance_, DType::Float);
		}

		if (dtype_ == DType::Fraction) {
			std::vector<std::vector<Rational>> out(rows_, std::vector<Rational>(1, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				Rational row_sum(0);
				for (std::size_t j = 0; j < cols_; ++j) {
					row_sum += frac_data_[i][j];
				}
				out[i][0] = row_sum;
			}
			return Matrix(out, tolerance_);
		}

		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(1, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[i][0] += float_data_[i][j];
			}
		}
		return Matrix(out, tolerance_, DType::Float);
	}

	Matrix concatenate(const Matrix& other, int index = 0) const {
		DType out_dtype = (dtype_ == DType::Fraction && other.dtype_ == DType::Fraction)
			? DType::Fraction
			: DType::Float;
		long double out_tol = std::max(tolerance_, other.tolerance_);

		if (index == 0) {
			if (cols_ != other.cols_) {
				throw std::invalid_argument("For row concatenation, column counts must match.");
			}

			if (out_dtype == DType::Fraction) {
				std::vector<std::vector<Rational>> out;
				out.reserve(rows_ + other.rows_);
				for (std::size_t i = 0; i < rows_; ++i) {
					out.push_back(frac_data_[i]);
				}
				for (std::size_t i = 0; i < other.rows_; ++i) {
					out.push_back(other.frac_data_[i]);
				}
				return Matrix(out, out_tol);
			}

			std::vector<std::vector<long double>> out;
			out.reserve(rows_ + other.rows_);
			for (std::size_t i = 0; i < rows_; ++i) {
				std::vector<long double> row(cols_, 0.0L);
				for (std::size_t j = 0; j < cols_; ++j) {
					row[j] = at_float(i, j);
				}
				out.push_back(row);
			}
			for (std::size_t i = 0; i < other.rows_; ++i) {
				std::vector<long double> row(cols_, 0.0L);
				for (std::size_t j = 0; j < cols_; ++j) {
					row[j] = other.at_float(i, j);
				}
				out.push_back(row);
			}
			return Matrix(out, out_tol, DType::Float);
		}

		if (index == 1) {
			if (rows_ != other.rows_) {
				throw std::invalid_argument("For column concatenation, row counts must match.");
			}

			if (out_dtype == DType::Fraction) {
				std::vector<std::vector<Rational>> out(rows_, std::vector<Rational>(cols_ + other.cols_, Rational(0)));
				for (std::size_t i = 0; i < rows_; ++i) {
					for (std::size_t j = 0; j < cols_; ++j) {
						out[i][j] = frac_data_[i][j];
					}
					for (std::size_t j = 0; j < other.cols_; ++j) {
						out[i][cols_ + j] = other.frac_data_[i][j];
					}
				}
				return Matrix(out, out_tol);
			}

			std::vector<std::vector<long double>> out(rows_, std::vector<long double>(cols_ + other.cols_, 0.0L));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					out[i][j] = at_float(i, j);
				}
				for (std::size_t j = 0; j < other.cols_; ++j) {
					out[i][cols_ + j] = other.at_float(i, j);
				}
			}
			return Matrix(out, out_tol, DType::Float);
		}

		throw std::invalid_argument("index must be 0 (rows) or 1 (cols).");
	}

	Matrix kronecker_product(const Matrix& other) const {
		DType out_dtype = (dtype_ == DType::Fraction && other.dtype_ == DType::Fraction)
			? DType::Fraction
			: DType::Float;
		long double out_tol = std::max(tolerance_, other.tolerance_);
		std::size_t out_rows = rows_ * other.rows_;
		std::size_t out_cols = cols_ * other.cols_;

		if (out_dtype == DType::Fraction) {
			std::vector<std::vector<Rational>> out(out_rows, std::vector<Rational>(out_cols, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					for (std::size_t p = 0; p < other.rows_; ++p) {
						for (std::size_t q = 0; q < other.cols_; ++q) {
							out[i * other.rows_ + p][j * other.cols_ + q] =
								frac_data_[i][j] * other.frac_data_[p][q];
						}
					}
				}
			}
			return Matrix(out, out_tol);
		}

		std::vector<std::vector<long double>> out(out_rows, std::vector<long double>(out_cols, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				for (std::size_t p = 0; p < other.rows_; ++p) {
					for (std::size_t q = 0; q < other.cols_; ++q) {
						out[i * other.rows_ + p][j * other.cols_ + q] =
							at_float(i, j) * other.at_float(p, q);
					}
				}
			}
		}
		return Matrix(out, out_tol, DType::Float);
	}

	Matrix Kronecker_product(const Matrix& other) const {
		return kronecker_product(other);
	}

	std::string str(int precision = 4) const {
		if (rows_ == 0) {
			return "[]";
		}
		if (cols_ == 0) {
			std::ostringstream oss;
			for (std::size_t i = 0; i < rows_; ++i) {
				oss << "[]";
				if (i + 1 < rows_) {
					oss << "\n";
				}
			}
			return oss.str();
		}

		std::ostringstream oss;
		oss << std::fixed << std::setprecision(precision);
		for (std::size_t i = 0; i < rows_; ++i) {
			oss << "[";
			for (std::size_t j = 0; j < cols_; ++j) {
				oss << std::setw(10) << static_cast<double>(at_float(i, j));
				if (j + 1 < cols_) {
					oss << " ";
				}
			}
			oss << "]";
			if (i + 1 < rows_) {
				oss << "\n";
			}
		}
		return oss.str();
	}

	void print(int precision = 4) const {
		std::cout << str(precision) << "\n";
	}

	friend std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
		os << mat.str();
		return os;
	}

	Matrix num_product(long long n) const {
		if (dtype_ == DType::Fraction) {
			Rational scalar(n);
			std::vector<std::vector<Rational>> out(rows_, std::vector<Rational>(cols_, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					out[i][j] = frac_data_[i][j] * scalar;
				}
			}
			return Matrix(out, tolerance_);
		}

		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(cols_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[i][j] = float_data_[i][j] * static_cast<long double>(n);
			}
		}
		return Matrix(out, tolerance_, DType::Float);
	}

	Matrix num_product(long double n) const {
		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(cols_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[i][j] = at_float(i, j) * n;
			}
		}
		return Matrix(out, tolerance_, DType::Float);
	}

	Matrix operator*(long long scalar) const {
		return num_product(scalar);
	}

	Matrix operator*(long double scalar) const {
		return num_product(scalar);
	}

	friend Matrix operator*(long long scalar, const Matrix& matrix) {
		return matrix.num_product(scalar);
	}

	friend Matrix operator*(long double scalar, const Matrix& matrix) {
		return matrix.num_product(scalar);
	}

	Matrix operator+(const Matrix& other) const {
		if (rows_ != other.rows_ || cols_ != other.cols_) {
			throw std::invalid_argument("Shapes must match for addition.");
		}

		DType out_dtype = (dtype_ == DType::Fraction && other.dtype_ == DType::Fraction)
			? DType::Fraction
			: DType::Float;
		long double out_tol = std::max(tolerance_, other.tolerance_);

		if (out_dtype == DType::Fraction) {
			std::vector<std::vector<Rational>> out(rows_, std::vector<Rational>(cols_, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					out[i][j] = frac_data_[i][j] + other.frac_data_[i][j];
				}
			}
			return Matrix(out, out_tol);
		}

		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(cols_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[i][j] = at_float(i, j) + other.at_float(i, j);
			}
		}
		return Matrix(out, out_tol, DType::Float);
	}

	Matrix operator-(const Matrix& other) const {
		if (rows_ != other.rows_ || cols_ != other.cols_) {
			throw std::invalid_argument("Shapes must match for subtraction.");
		}

		DType out_dtype = (dtype_ == DType::Fraction && other.dtype_ == DType::Fraction)
			? DType::Fraction
			: DType::Float;
		long double out_tol = std::max(tolerance_, other.tolerance_);

		if (out_dtype == DType::Fraction) {
			std::vector<std::vector<Rational>> out(rows_, std::vector<Rational>(cols_, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < cols_; ++j) {
					out[i][j] = frac_data_[i][j] - other.frac_data_[i][j];
				}
			}
			return Matrix(out, out_tol);
		}

		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(cols_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < cols_; ++j) {
				out[i][j] = at_float(i, j) - other.at_float(i, j);
			}
		}
		return Matrix(out, out_tol, DType::Float);
	}

	Matrix operator*(const Matrix& other) const {
		if (cols_ != other.rows_) {
			throw std::invalid_argument("Incompatible shapes for matrix multiplication.");
		}

		DType out_dtype = (dtype_ == DType::Fraction && other.dtype_ == DType::Fraction)
			? DType::Fraction
			: DType::Float;
		long double out_tol = std::max(tolerance_, other.tolerance_);

		if (out_dtype == DType::Fraction) {
			std::vector<std::vector<Rational>> out(rows_, std::vector<Rational>(other.cols_, Rational(0)));
			for (std::size_t i = 0; i < rows_; ++i) {
				for (std::size_t j = 0; j < other.cols_; ++j) {
					Rational val(0);
					for (std::size_t k = 0; k < cols_; ++k) {
						val += frac_data_[i][k] * other.frac_data_[k][j];
					}
					out[i][j] = val;
				}
			}
			return Matrix(out, out_tol);
		}

		std::vector<std::vector<long double>> out(rows_, std::vector<long double>(other.cols_, 0.0L));
		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < other.cols_; ++j) {
				long double val = 0.0L;
				for (std::size_t k = 0; k < cols_; ++k) {
					val += at_float(i, k) * other.at_float(k, j);
				}
				out[i][j] = val;
			}
		}
		return Matrix(out, out_tol, DType::Float);
	}

	Matrix power(int n) const {
		if (n < 0) {
			throw std::invalid_argument("Negative powers are not supported directly; use inverse() first.");
		}
		if (rows_ != cols_) {
			throw std::invalid_argument("Only square matrices can be exponentiated.");
		}

		Matrix result = Matrix::identity(rows_, tolerance_, dtype_);
		Matrix base = *this;
		int exp = n;
		while (exp > 0) {
			if (exp % 2 == 1) {
				result = result * base;
			}
			base = base * base;
			exp /= 2;
		}
		return result;
	}

	Matrix pow(int n) const {
		return power(n);
	}

	std::pair<Matrix, int> gauss_elimination(bool pivoting = true) const {
		Matrix mat = *this;
		int swap_count = 0;
		std::size_t r = 0;

		for (std::size_t c = 0; c < cols_ && r < rows_; ++c) {
			std::optional<std::size_t> pivot_row;

			if (pivoting) {
				long double best = 0.0L;
				for (std::size_t i = r; i < rows_; ++i) {
					if (mat.dtype_ == DType::Fraction) {
						if (!near_zero(mat.frac_data_[i][c], mat.tolerance_)) {
							long double score = std::fabs(mat.frac_data_[i][c].toLongDouble());
							if (score > best) {
								best = score;
								pivot_row = i;
							}
						}
					} else {
						if (!near_zero(mat.float_data_[i][c], mat.tolerance_)) {
							long double score = std::fabs(mat.float_data_[i][c]);
							if (score > best) {
								best = score;
								pivot_row = i;
							}
						}
					}
				}
			} else {
				for (std::size_t i = r; i < rows_; ++i) {
					if (mat.dtype_ == DType::Fraction) {
						if (!near_zero(mat.frac_data_[i][c], mat.tolerance_)) {
							pivot_row = i;
							break;
						}
					} else {
						if (!near_zero(mat.float_data_[i][c], mat.tolerance_)) {
							pivot_row = i;
							break;
						}
					}
				}
			}

			if (!pivot_row.has_value()) {
				continue;
			}

			if (pivot_row.value() != r) {
				mat.swap_rows(r, pivot_row.value());
				++swap_count;
			}

			if (mat.dtype_ == DType::Fraction) {
				Rational pivot = mat.frac_data_[r][c];
				for (std::size_t i = r + 1; i < rows_; ++i) {
					if (near_zero(mat.frac_data_[i][c], mat.tolerance_)) {
						continue;
					}
					Rational factor = mat.frac_data_[i][c] / pivot;
					for (std::size_t j = c; j < cols_; ++j) {
						mat.frac_data_[i][j] -= factor * mat.frac_data_[r][j];
						if (near_zero(mat.frac_data_[i][j], mat.tolerance_)) {
							mat.frac_data_[i][j] = Rational(0);
						}
					}
				}
			} else {
				long double pivot = mat.float_data_[r][c];
				for (std::size_t i = r + 1; i < rows_; ++i) {
					if (near_zero(mat.float_data_[i][c], mat.tolerance_)) {
						continue;
					}
					long double factor = mat.float_data_[i][c] / pivot;
					for (std::size_t j = c; j < cols_; ++j) {
						mat.float_data_[i][j] -= factor * mat.float_data_[r][j];
						if (near_zero(mat.float_data_[i][j], mat.tolerance_)) {
							mat.float_data_[i][j] = 0.0L;
						}
					}
				}
			}
			++r;
		}

		return {mat, swap_count};
	}

	long double determinant() const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Determinant is defined only for square matrices.");
		}

		auto reduced = gauss_elimination(true);
		const Matrix& ref = reduced.first;
		int swaps = reduced.second;

		if (dtype_ == DType::Fraction) {
			Rational det(1);
			for (std::size_t i = 0; i < rows_; ++i) {
				det *= ref.frac_data_[i][i];
			}
			if (swaps % 2 == 1) {
				det = -det;
			}
			if (near_zero(det, tolerance_)) {
				return 0.0L;
			}
			return det.toLongDouble();
		}

		long double det = 1.0L;
		for (std::size_t i = 0; i < rows_; ++i) {
			det *= ref.float_data_[i][i];
		}
		if (swaps % 2 == 1) {
			det = -det;
		}
		return near_zero(det, tolerance_) ? 0.0L : det;
	}

	long double det() const {
		return determinant();
	}

	Matrix inverse() const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Inverse is defined only for square matrices.");
		}

		std::size_t n = rows_;

		if (dtype_ == DType::Fraction) {
			std::vector<std::vector<Rational>> aug(n, std::vector<Rational>(2 * n, Rational(0)));
			for (std::size_t i = 0; i < n; ++i) {
				for (std::size_t j = 0; j < n; ++j) {
					aug[i][j] = frac_data_[i][j];
					aug[i][j + n] = (i == j) ? Rational(1) : Rational(0);
				}
			}

			for (std::size_t col = 0; col < n; ++col) {
				std::size_t pivot_row = col;
				long double best = std::fabs(aug[col][col].toLongDouble());
				for (std::size_t r = col + 1; r < n; ++r) {
					long double score = std::fabs(aug[r][col].toLongDouble());
					if (score > best) {
						best = score;
						pivot_row = r;
					}
				}

				if (near_zero(aug[pivot_row][col], tolerance_)) {
					throw std::runtime_error("Matrix is singular and cannot be inverted.");
				}

				if (pivot_row != col) {
					std::swap(aug[col], aug[pivot_row]);
				}

				Rational pivot = aug[col][col];
				for (std::size_t j = col; j < 2 * n; ++j) {
					aug[col][j] /= pivot;
				}

				for (std::size_t r = 0; r < n; ++r) {
					if (r == col) {
						continue;
					}
					Rational factor = aug[r][col];
					if (near_zero(factor, tolerance_)) {
						continue;
					}
					for (std::size_t j = col; j < 2 * n; ++j) {
						aug[r][j] -= factor * aug[col][j];
						if (near_zero(aug[r][j], tolerance_)) {
							aug[r][j] = Rational(0);
						}
					}
				}
			}

			std::vector<std::vector<Rational>> inv(n, std::vector<Rational>(n, Rational(0)));
			for (std::size_t i = 0; i < n; ++i) {
				for (std::size_t j = 0; j < n; ++j) {
					inv[i][j] = aug[i][j + n];
				}
			}
			return Matrix(inv, tolerance_);
		}

		std::vector<std::vector<long double>> aug(n, std::vector<long double>(2 * n, 0.0L));
		for (std::size_t i = 0; i < n; ++i) {
			for (std::size_t j = 0; j < n; ++j) {
				aug[i][j] = float_data_[i][j];
				aug[i][j + n] = (i == j) ? 1.0L : 0.0L;
			}
		}

		for (std::size_t col = 0; col < n; ++col) {
			std::size_t pivot_row = col;
			long double best = std::fabs(aug[col][col]);
			for (std::size_t r = col + 1; r < n; ++r) {
				long double score = std::fabs(aug[r][col]);
				if (score > best) {
					best = score;
					pivot_row = r;
				}
			}

			if (near_zero(aug[pivot_row][col], tolerance_)) {
				throw std::runtime_error("Matrix is singular and cannot be inverted.");
			}

			if (pivot_row != col) {
				std::swap(aug[col], aug[pivot_row]);
			}

			long double pivot = aug[col][col];
			for (std::size_t j = col; j < 2 * n; ++j) {
				aug[col][j] /= pivot;
			}

			for (std::size_t r = 0; r < n; ++r) {
				if (r == col) {
					continue;
				}
				long double factor = aug[r][col];
				if (near_zero(factor, tolerance_)) {
					continue;
				}
				for (std::size_t j = col; j < 2 * n; ++j) {
					aug[r][j] -= factor * aug[col][j];
					if (near_zero(aug[r][j], tolerance_)) {
						aug[r][j] = 0.0L;
					}
				}
			}
		}

		std::vector<std::vector<long double>> inv(n, std::vector<long double>(n, 0.0L));
		for (std::size_t i = 0; i < n; ++i) {
			for (std::size_t j = 0; j < n; ++j) {
				inv[i][j] = aug[i][j + n];
			}
		}
		return Matrix(inv, tolerance_, DType::Float);
	}

	int rank() const {
		auto reduced = gauss_elimination(true);
		const Matrix& ref = reduced.first;
		int rank_val = 0;
		for (std::size_t i = 0; i < ref.rows_; ++i) {
			bool has_nonzero = false;
			for (std::size_t j = 0; j < ref.cols_; ++j) {
				if (ref.dtype_ == DType::Fraction) {
					if (!near_zero(ref.frac_data_[i][j], ref.tolerance_)) {
						has_nonzero = true;
						break;
					}
				} else {
					if (!near_zero(ref.float_data_[i][j], ref.tolerance_)) {
						has_nonzero = true;
						break;
					}
				}
			}
			if (has_nonzero) {
				++rank_val;
			}
		}
		return rank_val;
	}

	std::tuple<Matrix, Matrix, Matrix> lu_decomposition(bool pivoting = true) const {
		if (rows_ != cols_) {
			throw std::invalid_argument("LU decomposition requires a square matrix.");
		}

		std::size_t n = rows_;
		if (dtype_ == DType::Fraction) {
			std::vector<std::vector<Rational>> U = frac_data_;
			std::vector<std::vector<Rational>> L(n, std::vector<Rational>(n, Rational(0)));
			std::vector<std::vector<Rational>> P(n, std::vector<Rational>(n, Rational(0)));

			for (std::size_t i = 0; i < n; ++i) {
				L[i][i] = Rational(1);
				P[i][i] = Rational(1);
			}

			for (std::size_t k = 0; k < n; ++k) {
				std::size_t pivot_row = k;
				if (pivoting) {
					long double best = std::fabs(U[k][k].toLongDouble());
					for (std::size_t r = k + 1; r < n; ++r) {
						long double score = std::fabs(U[r][k].toLongDouble());
						if (score > best) {
							best = score;
							pivot_row = r;
						}
					}
				}

				if (near_zero(U[pivot_row][k], tolerance_)) {
					throw std::runtime_error("Matrix is singular or LU requires pivoting at a zero pivot.");
				}

				if (pivot_row != k) {
					std::swap(U[k], U[pivot_row]);
					std::swap(P[k], P[pivot_row]);
					for (std::size_t c = 0; c < k; ++c) {
						std::swap(L[k][c], L[pivot_row][c]);
					}
				}

				for (std::size_t r = k + 1; r < n; ++r) {
					Rational factor = U[r][k] / U[k][k];
					L[r][k] = factor;
					for (std::size_t c = k; c < n; ++c) {
						U[r][c] -= factor * U[k][c];
						if (near_zero(U[r][c], tolerance_)) {
							U[r][c] = Rational(0);
						}
					}
				}
			}

			return {
				Matrix(P, tolerance_),
				Matrix(L, tolerance_),
				Matrix(U, tolerance_)
			};
		}

		std::vector<std::vector<long double>> U = float_data_;
		std::vector<std::vector<long double>> L(n, std::vector<long double>(n, 0.0L));
		std::vector<std::vector<long double>> P(n, std::vector<long double>(n, 0.0L));

		for (std::size_t i = 0; i < n; ++i) {
			L[i][i] = 1.0L;
			P[i][i] = 1.0L;
		}

		for (std::size_t k = 0; k < n; ++k) {
			std::size_t pivot_row = k;
			if (pivoting) {
				long double best = std::fabs(U[k][k]);
				for (std::size_t r = k + 1; r < n; ++r) {
					long double score = std::fabs(U[r][k]);
					if (score > best) {
						best = score;
						pivot_row = r;
					}
				}
			}

			if (near_zero(U[pivot_row][k], tolerance_)) {
				throw std::runtime_error("Matrix is singular or LU requires pivoting at a zero pivot.");
			}

			if (pivot_row != k) {
				std::swap(U[k], U[pivot_row]);
				std::swap(P[k], P[pivot_row]);
				for (std::size_t c = 0; c < k; ++c) {
					std::swap(L[k][c], L[pivot_row][c]);
				}
			}

			for (std::size_t r = k + 1; r < n; ++r) {
				long double factor = U[r][k] / U[k][k];
				L[r][k] = factor;
				for (std::size_t c = k; c < n; ++c) {
					U[r][c] -= factor * U[k][c];
					if (near_zero(U[r][c], tolerance_)) {
						U[r][c] = 0.0L;
					}
				}
			}
		}

		return {
			Matrix(P, tolerance_, DType::Float),
			Matrix(L, tolerance_, DType::Float),
			Matrix(U, tolerance_, DType::Float)
		};
	}

	std::pair<Matrix, Matrix> lu_decomposition_no_pivot() const {
		auto result = lu_decomposition(false);
		return {std::get<1>(result), std::get<2>(result)};
	}

	void luDecomposition(Matrix& L, Matrix& U) const {
		auto result = lu_decomposition(false);
		L = std::get<1>(result);
		U = std::get<2>(result);
	}

	std::pair<Matrix, Matrix> cholesky_decomposition() const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Cholesky decomposition requires a square matrix.");
		}
		if (!_is_symmetric()) {
			throw std::invalid_argument("Cholesky decomposition requires a symmetric matrix.");
		}

		std::size_t n = rows_;
		auto A = to_float_data();
		std::vector<std::vector<long double>> L(n, std::vector<long double>(n, 0.0L));

		for (std::size_t i = 0; i < n; ++i) {
			for (std::size_t j = 0; j <= i; ++j) {
				long double s = 0.0L;
				for (std::size_t k = 0; k < j; ++k) {
					s += L[i][k] * L[j][k];
				}

				if (i == j) {
					long double val = A[i][i] - s;
					if (val <= tolerance_) {
						throw std::runtime_error("Matrix is not positive definite.");
					}
					L[i][j] = std::sqrt(val);
				} else {
					if (std::fabs(L[j][j]) <= tolerance_) {
						throw std::runtime_error("Matrix is not positive definite.");
					}
					L[i][j] = (A[i][j] - s) / L[j][j];
				}
			}
		}

		Matrix L_mat(L, tolerance_, DType::Float);
		return {L_mat, L_mat.transpose()};
	}

	Matrix gram_schmidt(const std::string& method = "mgs") const {
		if (rows_ == 0 || cols_ == 0) {
			throw std::invalid_argument("Cannot run Gram-Schmidt on an empty matrix.");
		}

		std::string method_key = to_lower(method);
		if (
			method_key != "cgs" && method_key != "classical" &&
			method_key != "mgs" && method_key != "modified"
		) {
			throw std::invalid_argument("method must be 'cgs'/'classical' or 'mgs'/'modified'.");
		}

		auto A = to_float_data();
		std::vector<std::vector<long double>> cols(cols_, std::vector<long double>(rows_, 0.0L));
		for (std::size_t j = 0; j < cols_; ++j) {
			for (std::size_t i = 0; i < rows_; ++i) {
				cols[j][i] = A[i][j];
			}
		}

		std::vector<std::vector<long double>> q_cols;

		if (method_key == "cgs" || method_key == "classical") {
			for (const auto& col : cols) {
				std::vector<long double> v = col;
				for (const auto& q : q_cols) {
					long double proj = dot(col, q);
					for (std::size_t i = 0; i < rows_; ++i) {
						v[i] -= proj * q[i];
					}
				}
				long double nrm = norm(v);
				if (nrm <= tolerance_) {
					throw std::runtime_error("Columns are linearly dependent; CGS failed.");
				}
				for (long double& x : v) {
					x /= nrm;
				}
				q_cols.push_back(v);
			}
		} else {
			for (const auto& col : cols) {
				std::vector<long double> v = col;
				for (const auto& q : q_cols) {
					long double proj = dot(v, q);
					for (std::size_t i = 0; i < rows_; ++i) {
						v[i] -= proj * q[i];
					}
				}
				long double nrm = norm(v);
				if (nrm <= tolerance_) {
					throw std::runtime_error("Columns are linearly dependent; MGS failed.");
				}
				for (long double& x : v) {
					x /= nrm;
				}
				q_cols.push_back(v);
			}
		}

		std::vector<std::vector<long double>> q_data(rows_, std::vector<long double>(cols_, 0.0L));
		for (std::size_t j = 0; j < cols_; ++j) {
			for (std::size_t i = 0; i < rows_; ++i) {
				q_data[i][j] = q_cols[j][i];
			}
		}

		return Matrix(q_data, tolerance_, DType::Float);
	}

	std::pair<Matrix, Matrix> _qr_householder() const {
		if (rows_ == 0 || cols_ == 0) {
			throw std::invalid_argument("Cannot run QR on an empty matrix.");
		}

		auto R = to_float_data();
		auto Q = Matrix::identity(rows_, tolerance_, DType::Float).to_float_data();
		std::size_t k_max = std::min(rows_, cols_);

		for (std::size_t k = 0; k < k_max; ++k) {
			std::vector<long double> x(rows_ - k, 0.0L);
			for (std::size_t i = k; i < rows_; ++i) {
				x[i - k] = R[i][k];
			}
			long double norm_x = norm(x);
			if (norm_x <= tolerance_) {
				continue;
			}

			long double sign = (x[0] < 0.0L) ? -1.0L : 1.0L;
			std::vector<long double> v = x;
			v[0] += sign * norm_x;
			long double norm_v = norm(v);
			if (norm_v <= tolerance_) {
				continue;
			}
			for (long double& val : v) {
				val /= norm_v;
			}

			for (std::size_t j = k; j < cols_; ++j) {
				long double dot_vr = 0.0L;
				for (std::size_t i = 0; i < v.size(); ++i) {
					dot_vr += v[i] * R[k + i][j];
				}
				for (std::size_t i = 0; i < v.size(); ++i) {
					R[k + i][j] -= 2.0L * v[i] * dot_vr;
				}
			}

			for (std::size_t i = 0; i < rows_; ++i) {
				long double dot_qv = 0.0L;
				for (std::size_t t = 0; t < v.size(); ++t) {
					dot_qv += v[t] * Q[i][k + t];
				}
				for (std::size_t t = 0; t < v.size(); ++t) {
					Q[i][k + t] -= 2.0L * v[t] * dot_qv;
				}
			}
		}

		for (std::size_t i = 0; i < rows_; ++i) {
			for (std::size_t j = 0; j < std::min(i, cols_); ++j) {
				if (std::fabs(R[i][j]) <= tolerance_) {
					R[i][j] = 0.0L;
				}
			}
		}

		return {
			Matrix(Q, tolerance_, DType::Float),
			Matrix(R, tolerance_, DType::Float)
		};
	}

	std::pair<Matrix, Matrix> qr_decomposition(const std::string& method = "householder") const {
		std::string method_key = to_lower(method);

		if (method_key == "householder" || method_key == "hh") {
			return _qr_householder();
		}

		if (
			method_key == "mgs" || method_key == "modified" ||
			method_key == "cgs" || method_key == "classical"
		) {
			Matrix Q = gram_schmidt(method_key);
			Matrix A_float(to_float_data(), tolerance_, DType::Float);
			Matrix R = Q.transpose() * A_float;

			auto r_data = R.to_float_data();
			for (std::size_t i = 0; i < R.rows_; ++i) {
				for (std::size_t j = 0; j < R.cols_; ++j) {
					if (i > j && std::fabs(r_data[i][j]) <= tolerance_) {
						r_data[i][j] = 0.0L;
					}
				}
			}

			return {
				Q,
				Matrix(r_data, tolerance_, DType::Float)
			};
		}

		throw std::invalid_argument("Unknown QR method.");
	}

	std::pair<Matrix, Matrix> QR_decomposition(const std::string& method = "householder") const {
		return qr_decomposition(method);
	}

	void qrDecomposition(Matrix& Q, Matrix& R) const {
		auto qr = qr_decomposition("householder");
		Q = qr.first;
		R = qr.second;
	}

	std::vector<long double> eigenvalues_qr(
		int iterations = 100,
		const std::string& qr_method = "householder"
	) const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Eigenvalues are defined only for square matrices.");
		}

		Matrix current(to_float_data(), tolerance_, DType::Float);
		for (int i = 0; i < iterations; ++i) {
			auto qr = current.qr_decomposition(qr_method);
			current = qr.second * qr.first;
		}

		std::vector<long double> vals(rows_, 0.0L);
		for (std::size_t i = 0; i < rows_; ++i) {
			long double v = current.at_float(i, i);
			vals[i] = std::fabs(v) <= tolerance_ ? 0.0L : v;
		}
		return vals;
	}

	std::pair<std::vector<long double>, Matrix> eigenpairs_jacobi(
		int max_iterations = 1000,
		bool sort = true
	) const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Eigenvalues are defined only for square matrices.");
		}
		if (!_is_symmetric()) {
			throw std::invalid_argument("Jacobi method requires a symmetric matrix.");
		}

		std::size_t n = rows_;
		auto A = to_float_data();
		auto V = Matrix::identity(n, tolerance_, DType::Float).to_float_data();

		for (int iter = 0; iter < max_iterations; ++iter) {
			long double max_val = 0.0L;
			std::size_t p = 0;
			std::size_t q = 0;
			for (std::size_t i = 0; i < n; ++i) {
				for (std::size_t j = i + 1; j < n; ++j) {
					long double candidate = std::fabs(A[i][j]);
					if (candidate > max_val) {
						max_val = candidate;
						p = i;
						q = j;
					}
				}
			}

			if (max_val < tolerance_) {
				break;
			}

			long double tau = (A[q][q] - A[p][p]) / (2.0L * A[p][q]);
			long double sign_tau = (tau >= 0.0L) ? 1.0L : -1.0L;
			long double t = sign_tau / (std::fabs(tau) + std::sqrt(tau * tau + 1.0L));
			long double c = 1.0L / std::sqrt(1.0L + t * t);
			long double s = t * c;

			for (std::size_t i = 0; i < n; ++i) {
				if (i != p && i != q) {
					long double aip = A[i][p];
					long double aiq = A[i][q];
					A[i][p] = c * aip - s * aiq;
					A[p][i] = A[i][p];
					A[i][q] = c * aiq + s * aip;
					A[q][i] = A[i][q];
				}
			}

			long double app = A[p][p];
			long double aqq = A[q][q];
			long double apq = A[p][q];
			A[p][p] = c * c * app - 2.0L * s * c * apq + s * s * aqq;
			A[q][q] = s * s * app + 2.0L * s * c * apq + c * c * aqq;
			A[p][q] = 0.0L;
			A[q][p] = 0.0L;

			for (std::size_t i = 0; i < n; ++i) {
				long double vip = V[i][p];
				long double viq = V[i][q];
				V[i][p] = c * vip - s * viq;
				V[i][q] = c * viq + s * vip;
			}
		}

		std::vector<long double> eigenvalues(n, 0.0L);
		for (std::size_t i = 0; i < n; ++i) {
			eigenvalues[i] = std::fabs(A[i][i]) > tolerance_ ? A[i][i] : 0.0L;
		}

		if (sort) {
			std::vector<std::size_t> idx(n, 0);
			std::iota(idx.begin(), idx.end(), 0);
			std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
				return eigenvalues[a] > eigenvalues[b];
			});

			std::vector<long double> sorted_vals(n, 0.0L);
			std::vector<std::vector<long double>> sorted_V(n, std::vector<long double>(n, 0.0L));
			for (std::size_t i = 0; i < n; ++i) {
				sorted_vals[i] = eigenvalues[idx[i]];
				for (std::size_t r = 0; r < n; ++r) {
					sorted_V[r][i] = V[r][idx[i]];
				}
			}
			eigenvalues = sorted_vals;
			V = sorted_V;
		}

		return {eigenvalues, Matrix(V, tolerance_, DType::Float)};
	}

	std::pair<std::vector<long double>, Matrix> eigenpairs(
		const std::string& method = "auto",
		int iterations = 100,
		int max_iterations = 1000,
		const std::string& qr_method = "householder"
	) const {
		std::string method_key = to_lower(method);
		if (method_key == "auto") {
			method_key = _is_symmetric() ? "jacobi" : "qr";
		}

		if (method_key == "jacobi") {
			return eigenpairs_jacobi(max_iterations, true);
		}

		if (method_key == "qr") {
			return {
				eigenvalues_qr(iterations, qr_method),
				Matrix()
			};
		}

		throw std::invalid_argument("method must be 'auto', 'jacobi', or 'qr'.");
	}

	std::vector<long double> eigenvalues(
		const std::string& method = "auto",
		int iterations = 100,
		int max_iterations = 1000,
		const std::string& qr_method = "householder"
	) const {
		return eigenpairs(method, iterations, max_iterations, qr_method).first;
	}

	std::tuple<Matrix, Matrix, Matrix> singular_value_decomposition() const {
		std::size_t m = rows_;
		std::size_t n = cols_;

		auto A_float = to_float_data();
		Matrix AtA((this->transpose() * (*this)).to_float_data(), tolerance_, DType::Float);
		auto eig = AtA.eigenpairs_jacobi(2000, true);
		std::vector<long double> eigvals = eig.first;
		Matrix V = eig.second;

		for (long double& ev : eigvals) {
			ev = ev > tolerance_ ? ev : 0.0L;
		}

		std::vector<long double> singular_values(eigvals.size(), 0.0L);
		for (std::size_t i = 0; i < eigvals.size(); ++i) {
			singular_values[i] = std::sqrt(eigvals[i]);
		}

		std::size_t k = std::min(m, n);
		std::vector<std::vector<long double>> sigma_data(m, std::vector<long double>(n, 0.0L));
		for (std::size_t i = 0; i < k; ++i) {
			sigma_data[i][i] = singular_values[i];
		}

		auto V_data = V.to_float_data();
		std::vector<std::vector<long double>> v_cols(n, std::vector<long double>(n, 0.0L));
		for (std::size_t c = 0; c < n; ++c) {
			for (std::size_t r = 0; r < n; ++r) {
				v_cols[c][r] = V_data[r][c];
			}
		}

		std::vector<std::vector<long double>> u_cols(k);
		for (std::size_t i = 0; i < k; ++i) {
			long double sigma = singular_values[i];
			if (sigma <= tolerance_) {
				continue;
			}
			auto Av = mat_vec_mul(A_float, v_cols[i]);
			std::vector<long double> candidate(Av.size(), 0.0L);
			for (std::size_t t = 0; t < Av.size(); ++t) {
				candidate[t] = Av[t] / sigma;
			}
			u_cols[i] = normalize_vec(candidate, tolerance_);
		}

		std::vector<std::vector<long double>> known_u;
		for (const auto& col : u_cols) {
			if (!col.empty()) {
				known_u.push_back(col);
			}
		}
		auto full_basis = complete_orthonormal_basis(known_u, m, tolerance_);

		std::size_t basis_cursor = known_u.size();
		for (std::size_t i = 0; i < k; ++i) {
			if (u_cols[i].empty()) {
				u_cols[i] = full_basis[basis_cursor++];
			}
		}

		while (u_cols.size() < m) {
			u_cols.push_back(full_basis[basis_cursor++]);
		}

		std::vector<std::vector<long double>> U_data(m, std::vector<long double>(m, 0.0L));
		for (std::size_t row = 0; row < m; ++row) {
			for (std::size_t col = 0; col < m; ++col) {
				U_data[row][col] = u_cols[col][row];
			}
		}

		Matrix U(U_data, tolerance_, DType::Float);
		Matrix Sigma(sigma_data, tolerance_, DType::Float);
		Matrix V_T(V_data, tolerance_, DType::Float);
		V_T = V_T.transpose();
		return {U, Sigma, V_T};
	}

	std::tuple<Matrix, Matrix, Matrix> svd() const {
		return singular_value_decomposition();
	}

	Matrix jordan_normal_form() const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Jordan normal form is defined only for square matrices.");
		}
		if (rows_ == 0) {
			return Matrix();
		}

		std::size_t n = rows_;
		Matrix A_mat = to_fraction_matrix();

		std::vector<Rational> c(n + 1, Rational(0));
		c[n] = Rational(1);

		Matrix M_mat(
			std::vector<std::vector<Rational>>(n, std::vector<Rational>(n, Rational(0))),
			tolerance_
		);

		for (std::size_t k = 1; k <= n; ++k) {
			Rational c_prev = c[n - k + 1];

			std::vector<std::vector<Rational>> c_prev_I(n, std::vector<Rational>(n, Rational(0)));
			for (std::size_t i = 0; i < n; ++i) {
				c_prev_I[i][i] = c_prev;
			}

			M_mat = (A_mat * M_mat) + Matrix(c_prev_I, tolerance_);
			Matrix AM = A_mat * M_mat;

			Rational trace_AM(0);
			for (std::size_t i = 0; i < n; ++i) {
				trace_AM += AM.at_fraction(i, i);
			}

			c[n - k] = Rational(-1, static_cast<std::int64_t>(k)) * trace_AM;
		}

		std::map<Rational, int> roots;
		std::vector<Rational> poly = c;

		while (poly.size() > 1 && poly.front().isZero()) {
			roots[Rational(0)] += 1;
			poly.erase(poly.begin());
		}

		if (poly.size() > 1) {
			auto lcm64 = [](std::int64_t a, std::int64_t b) -> std::int64_t {
				a = std::llabs(a);
				b = std::llabs(b);
				if (a == 0 || b == 0) {
					return 0;
				}
				return (a / std::gcd(a, b)) * b;
			};

			std::int64_t lcm_den = 1;
			for (const auto& coeff : poly) {
				lcm_den = lcm64(lcm_den, coeff.denominator());
			}

			std::vector<Rational> int_coeffs;
			int_coeffs.reserve(poly.size());
			for (const auto& coeff : poly) {
				int_coeffs.push_back(
					Rational(coeff.numerator() * (lcm_den / coeff.denominator()), 1)
				);
			}

			std::int64_t c0 = int_coeffs.front().numerator();
			std::int64_t cn = int_coeffs.back().numerator();

			auto p_candidates = get_factors(c0);
			auto q_candidates = get_factors(cn);

			std::map<Rational, bool> candidates;
			for (std::int64_t p : p_candidates) {
				for (std::int64_t q : q_candidates) {
					candidates[Rational(p, q)] = true;
					candidates[Rational(-p, q)] = true;
				}
			}

			for (const auto& entry : candidates) {
				const Rational& cand = entry.first;
				while (int_coeffs.size() > 1 && evaluate_poly(int_coeffs, cand).isZero()) {
					roots[cand] += 1;
					int_coeffs = synthetic_div(int_coeffs, cand);
				}
			}

			if (int_coeffs.size() > 1) {
				throw std::runtime_error("Matrix has non-rational eigenvalues; exact Jordan form is unsupported here.");
			}
		}

		std::vector<std::pair<Rational, int>> jordan_blocks;
		for (const auto& root : roots) {
			const Rational& eigenvalue = root.first;
			int alg_mult = root.second;

			std::vector<std::vector<Rational>> diag_B_data(n, std::vector<Rational>(n, Rational(0)));
			for (std::size_t i = 0; i < n; ++i) {
				diag_B_data[i][i] = -eigenvalue;
			}
			Matrix B = A_mat + Matrix(diag_B_data, tolerance_);

			std::vector<int> ranks;
			ranks.push_back(static_cast<int>(n));
			Matrix current_B = Matrix::identity(n, tolerance_, DType::Fraction);

			for (int k = 1; k <= alg_mult + 1; ++k) {
				current_B = current_B * B;
				ranks.push_back(current_B.rank());

				if (ranks.size() >= 3 && ranks[ranks.size() - 1] == ranks[ranks.size() - 2]) {
					while (static_cast<int>(ranks.size()) <= alg_mult + 2) {
						ranks.push_back(ranks.back());
					}
					break;
				}
			}

			for (int k = 1; k <= alg_mult; ++k) {
				if (k + 1 < static_cast<int>(ranks.size())) {
					int num_blocks = ranks[k - 1] - 2 * ranks[k] + ranks[k + 1];
					for (int b = 0; b < num_blocks; ++b) {
						jordan_blocks.push_back({eigenvalue, k});
					}
				}
			}
		}

		std::vector<std::vector<Rational>> J_data(n, std::vector<Rational>(n, Rational(0)));
		std::size_t idx = 0;
		for (const auto& block : jordan_blocks) {
			const Rational& eigenvalue = block.first;
			int block_size = block.second;
			for (int i = 0; i < block_size && idx + static_cast<std::size_t>(i) < n; ++i) {
				J_data[idx + static_cast<std::size_t>(i)][idx + static_cast<std::size_t>(i)] = eigenvalue;
				if (i < block_size - 1 && idx + static_cast<std::size_t>(i + 1) < n) {
					J_data[idx + static_cast<std::size_t>(i)][idx + static_cast<std::size_t>(i + 1)] = Rational(1);
				}
			}
			idx += static_cast<std::size_t>(block_size);
			if (idx >= n) {
				break;
			}
		}

		return Matrix(J_data, tolerance_);
	}

	Matrix jordanNormalForm() const {
		return jordan_normal_form();
	}

	void dominantEigen(long double& eigenvalue, Matrix& eigenvector, int iterations = 1000) const {
		if (rows_ != cols_) {
			throw std::invalid_argument("Matrix must be square.");
		}

		eigenvector = Matrix(rows_, 1, 1.0L, tolerance_, DType::Float);
		for (int it = 0; it < iterations; ++it) {
			Matrix next_vec = (*this) * eigenvector;
			long double nrm = 0.0L;
			for (std::size_t i = 0; i < rows_; ++i) {
				long double val = next_vec.at_float(i, 0);
				nrm += val * val;
			}
			nrm = std::sqrt(nrm);
			if (nrm <= tolerance_) {
				throw std::runtime_error("Power iteration failed due to near-zero norm.");
			}
			for (std::size_t i = 0; i < rows_; ++i) {
				eigenvector.set_value(i, 0, next_vec.at_float(i, 0) / nrm);
			}
		}

		Matrix Ax = (*this) * eigenvector;
		eigenvalue = 0.0L;
		for (std::size_t i = 0; i < rows_; ++i) {
			eigenvalue += eigenvector.at_float(i, 0) * Ax.at_float(i, 0);
		}
	}

	void fullRankDecomposition(Matrix& F, Matrix& G) const {
		Matrix rref(to_float_data(), tolerance_, DType::Float);
		std::vector<std::size_t> pivot_cols;
		std::size_t lead = 0;

		for (std::size_t r = 0; r < rows_; ++r) {
			if (lead >= cols_) {
				break;
			}

			std::size_t i = r;
			while (i < rows_ && std::fabs(rref.at_float(i, lead)) <= tolerance_) {
				++i;
			}
			if (i == rows_) {
				++lead;
				--r;
				continue;
			}

			pivot_cols.push_back(lead);
			rref.swap_rows(i, r);

			long double lv = rref.at_float(r, lead);
			for (std::size_t c = 0; c < cols_; ++c) {
				rref.set_value(r, c, rref.at_float(r, c) / lv);
			}

			for (std::size_t i2 = 0; i2 < rows_; ++i2) {
				if (i2 == r) {
					continue;
				}
				long double lv2 = rref.at_float(i2, lead);
				for (std::size_t c = 0; c < cols_; ++c) {
					rref.set_value(i2, c, rref.at_float(i2, c) - lv2 * rref.at_float(r, c));
				}
			}

			++lead;
		}

		std::size_t rank_val = pivot_cols.size();
		F = Matrix(rows_, rank_val, 0.0L, tolerance_, DType::Float);
		G = Matrix(rank_val, cols_, 0.0L, tolerance_, DType::Float);

		for (std::size_t i = 0; i < rank_val; ++i) {
			for (std::size_t r = 0; r < rows_; ++r) {
				F.set_value(r, i, at_float(r, pivot_cols[i]));
			}
			for (std::size_t c = 0; c < cols_; ++c) {
				G.set_value(i, c, rref.at_float(i, c));
			}
		}
	}
};

Matrix I(std::size_t n, long double tolerance = 1e-10L, const std::string& dtype = "fraction") {
	return Matrix::identity(n, tolerance, parse_dtype(dtype));
}

Matrix concatenate(const std::vector<Matrix>& items, int axis = 0) {
	if (items.empty()) {
		throw std::invalid_argument("No matrices to concatenate.");
	}

	Matrix result = items[0].copy();
	for (std::size_t i = 1; i < items.size(); ++i) {
		result = result.concatenate(items[i], axis);
	}
	return result;
}

#ifdef MATLIBRARY_DEMO
int main() {
	std::cout << "--- Matrix Library Demo ---\n\n";

	Matrix A(3, 3, 0.0L, 1e-10L, DType::Float);
	A(0, 0) = -1.0L; A(0, 1) = 1.0L; A(0, 2) = 0.0L;
	A(1, 0) = -4.0L; A(1, 1) = 3.0L; A(1, 2) = 0.0L;
	A(2, 0) = 1.0L;  A(2, 1) = 0.0L; A(2, 2) = 2.0L;

	std::cout << "Matrix A:\n";
	A.print();

	std::cout << "Determinant of A: " << static_cast<double>(A.determinant()) << "\n\n";

	std::cout << "QR Eigenvalue Approximation:\n";
	auto evals = A.eigenvalues_qr(120);
	for (std::size_t i = 0; i < evals.size(); ++i) {
		std::cout << "lambda[" << i << "] = " << static_cast<double>(evals[i]) << "\n";
	}

	std::cout << "\n--- End ---\n";
	return 0;
}
#endif