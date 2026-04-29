#pragma once

#ifndef __MATLIB__
#include "matlib.hh"
#endif

#include <algorithm>
#include <execution>
#include <format>
#include <stdexcept>

#include "config.hh"

namespace MATLIB {
template <typename Scalar>
    requires concepts::is_ring<Scalar>
class Matrix<Scalar>::check {
   public:
    static bool shape_same(Shape shape_1, Shape shape_2) {
        if constexpr (config::kStrict) {
            if (shape_1 != shape_2) {
                throw std::runtime_error(std::format(
                    "Shape does not match! ({}, {}) != ({}, {})", shape_1.rows,
                    shape_1.cols, shape_2.rows, shape_2.cols));
                return false;
            } else {
                return true;
            }
        } else {
            return shape_1 == shape_2;
        }
    }

    static bool shape_suit(Shape shape_1, Shape shape_2) {
        if constexpr (config::kStrict) {
            if (shape_1.cols != shape_2.rows) {
                throw std::runtime_error(std::format(
                    "Shape does not match! ({}, {}) * ({}, {}) fails!",
                    shape_1.rows, shape_1.cols, shape_2.rows, shape_2.cols));
                return false;
            } else {
                return true;
            }
        } else {
            return shape_1.cols == shape_2.rows;
        }
    }
};

template <typename Scalar>
    requires concepts::is_ring<Scalar>
Matrix<Scalar>::Matrix(Scalar value, ::std::size_t rows, ::std::size_t cols)
    : shape_(rows, cols), data_(rows * cols, value) {}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
Matrix<Scalar>::Matrix(Scalar value, Shape shape)
    : shape_(shape), data_(shape.rows * shape.cols, value) {}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
[[nodiscard]] Scalar& Matrix<Scalar>::operator()(::std::size_t row,
                                                 ::std::size_t col) {
    return data_[row * shape_.cols + col];
}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
[[nodiscard]] const Scalar Matrix<Scalar>::operator()(std::size_t row,
                                                      std::size_t col) const {
    return data_[row * shape_.cols + col];
}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
void Matrix<Scalar>::reshape(::std::size_t new_rows, ::std::size_t new_cols) {
    if (new_rows * new_cols != shape_.rows * shape_.cols) {
        throw std::runtime_error(
            std::format("Shape does not match! {} * {} != {} * {}", new_rows,
                        new_cols, shape_.cols, shape_.rows));
    }
    shape_ = Shape(new_rows, new_cols);
}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
Matrix<Scalar>::Shape Matrix<Scalar>::shape() const {
    return shape_;
}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
Matrix<Scalar> operator+(const Matrix<Scalar>& A, const Matrix<Scalar>& B) {
    Matrix<Scalar>::check::shape_same(A.shape(), B.shape());
    Matrix<Scalar> C(A.shape());

    if constexpr (config::kAlwaysParallel) {
        ::std::transform(::std::execution::par, A.data_.begin(), A.data_.end(),
                         B.data_.begin(), C.data_.begin(),
                         [](Scalar a, Scalar b) { return a + b; });
    } else {
        if (C.data_.size() <= config::kMinParallel) {
            for (::std::size_t i = 0; i < C.data_.size(); i++) {
                C.data_[i] = A.data_[i] + B.data_[i];
            }
        } else {
            ::std::transform(::std::execution::par, A.data_.begin(),
                             A.data_.end(), B.data_.begin(), C.data_.begin(),
                             [](Scalar a, Scalar b) { return a + b; });
        }
    }
    return C;
}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
Matrix<Scalar> operator-(const Matrix<Scalar>& A, const Matrix<Scalar>& B) {
    Matrix<Scalar>::check::shape_same(A.shape(), B.shape());
    Matrix<Scalar> C(A.shape());

    if constexpr (config::kAlwaysParallel) {
        ::std::transform(::std::execution::par, A.data_.begin(), A.data_.end(),
                         B.data_.begin(), C.data_.begin(),
                         [](Scalar a, Scalar b) { return a - b; });
    } else {
        if (C.data_.size() <= config::kMinParallel) {
            for (::std::size_t i = 0; i < C.data_.size(); i++) {
                C.data_[i] = A.data_[i] - B.data_[i];
            }
        } else {
            ::std::transform(::std::execution::par, A.data_.begin(),
                             A.data_.end(), B.data_.begin(), C.data_.begin(),
                             [](Scalar a, Scalar b) { return a - b; });
        }
    }
    return C;
}

template <typename Scalar>
    requires concepts::is_ring<Scalar>
Matrix<Scalar> operator*(const Matrix<Scalar>& A, const Matrix<Scalar>& B) {
    Matrix<Scalar>::check::shape_suit(A.shape(), B.shape());

    const auto M = A.shape().rows, K = A.shape().cols, N = B.shape().cols;

    Matrix<Scalar> C(Scalar{0}, {A.shape().rows, B.shape().cols});

    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t t = 0; t < K; ++t) {
            const Scalar a = A.data_[i * K + t];

            for (std::size_t j = 0; j < N; ++j) {
                C.data_[i * N + j] =
                    C.data_[i * N + j] + a * B.data_[t * N + j];
            }
        }
    }

    return C;
}

template <typename Scalar>
bool operator==(const Matrix<Scalar>& A, const Matrix<Scalar>& B) {
    Matrix<Scalar>::check::shape_same(A.shape(), B.shape());

    for (std::size_t i = 0; i < A.data_.size(); i++) {
        if (A.data_[i] != B.data_[i]) [[unlikely]] {
            return false;
        }
    }

    return true;
}
}  // namespace MATLIB