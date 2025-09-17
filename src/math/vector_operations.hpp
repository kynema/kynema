#pragma once

#include <Eigen/Dense>
#include <Kokkos_Core.hpp>

namespace kynema::math {

/// Converts a 3x1 vector to a 3x3 skew-symmetric matrix and returns the result
template <typename VectorType, typename MatrixType>
KOKKOS_INLINE_FUNCTION void VecTilde(const VectorType& vector, const MatrixType& matrix) {
    matrix(0, 0) = 0.;
    matrix(0, 1) = -vector(2);
    matrix(0, 2) = vector(1);
    matrix(1, 0) = vector(2);
    matrix(1, 1) = 0.;
    matrix(1, 2) = -vector(0);
    matrix(2, 0) = -vector(1);
    matrix(2, 1) = vector(0);
    matrix(2, 2) = 0.;
}

/// Calculate the dot product between two vector views
template <typename AVectorType, typename BVectorType>
KOKKOS_INLINE_FUNCTION double DotProduct(const AVectorType& a, const BVectorType& b) {
    double sum = 0.;
    for (auto i = 0; i < a.extent_int(0); ++i) {
        sum += a(i) * b(i);
    }
    return sum;
}

/// Calculate the cross product between two vector views
template <typename VectorType>
KOKKOS_INLINE_FUNCTION void CrossProduct(
    const VectorType& a, const VectorType& b, const VectorType& c
) {
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
}

}  // namespace kynema::math
