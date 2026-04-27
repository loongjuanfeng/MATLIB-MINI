#pragma once

#include <cstddef>
#include <vector>

#include "concepts.hh"

namespace MATLIB {
template <typename Scalar>
    requires concepts::is_ring<Scalar>
class Matrix {
   public:
    struct Shape {
        std::size_t rows;
        std::size_t cols;

        Shape() = delete;
        Shape(std::size_t rows_, std::size_t cols_)
            : rows(rows_), cols(cols_) {}
        friend bool operator==(Shape shape_1, Shape shape_2) {
            return shape_1.rows == shape_2.rows && shape_1.cols == shape_2.cols;
        }
    };

   private:
    class check;

   private:
    Shape shape_;
    std::vector<Scalar> data_;
    Matrix(Shape shape) : shape_(shape), data_(shape.rows * shape.cols) {}

   public:
    Matrix() = delete;
    Matrix(Scalar value, std::size_t rows, std::size_t cols);
    Matrix(Scalar value, Shape shape);

    [[nodiscard]] const Scalar operator()(std::size_t row,
                                          std::size_t col) const;
    [[nodiscard]] Scalar& operator()(::std::size_t row, ::std::size_t col);

    void reshape(std::size_t new_rows, std::size_t new_cols);
    Shape shape() const;

    template <typename T>
        requires concepts::is_ring<T>
    friend Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B);

    template <typename T>
        requires concepts::is_ring<T>
    friend Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B);

    template <typename T>
        requires concepts::is_ring<T>
    friend Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B);

    template <typename T>
    friend bool operator==(const Matrix<T>& A, const Matrix<T>& B);
};
}  // namespace MATLIB

#include "matlib.cc"
