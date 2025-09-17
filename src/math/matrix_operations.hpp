#pragma once

#include <array>

#include <Eigen/Geometry>
#include <Kokkos_Core.hpp>

namespace kynema::math {

/**
 * @brief Computes AX(A) of a square matrix
 * @details AX(A) = tr(A)/2 * I - A/2, where I is the identity matrix
 * @param A Input square matrix
 * @param AX_A Output matrix containing the result
 */
template <typename Matrix>
KOKKOS_INLINE_FUNCTION void AX_Matrix(const Matrix& A, const Matrix& AX_A) {
    double trace{0.};
    for (auto i = 0; i < A.extent_int(0); ++i) {
        trace += A(i, i);
    }

    trace /= 2.;

    for (auto i = 0; i < A.extent_int(0); ++i) {
        for (auto j = 0; j < A.extent_int(1); ++j) {
            AX_A(i, j) = -A(i, j) / 2.;
        }
        AX_A(i, i) += trace;
    }
}

/**
 * @brief Computes the axial vector (also known as the vector representation) of a 3x3 skew-symmetric
 * matrix
 * @details The axial vector is defined as [w₁, w₂, w₃]ᵀ where:
 *          w₁ = (m₃₂ - m₂₃)/2
 *          w₂ = (m₁₃ - m₃₁)/2
 *          w₃ = (m₂₁ - m₁₂)/2
 * @param m Input 3x3 rotation matrix
 * @param v Output vector to store the result
 * @pre Matrix m must be 3x3 and vector v must have size 3
 */
template <typename Matrix, typename Vector>
KOKKOS_INLINE_FUNCTION void AxialVectorOfMatrix(const Matrix& m, const Vector& v) {
    v(0) = (m(2, 1) - m(1, 2)) / 2.;
    v(1) = (m(0, 2) - m(2, 0)) / 2.;
    v(2) = (m(1, 0) - m(0, 1)) / 2.;
}

inline std::array<std::array<double, 6>, 6> RotateMatrix6(
    const std::array<std::array<double, 6>, 6>& m, const std::array<double, 4>& q
) {
    const auto quat = Eigen::Quaternion<double>(q[0], q[1], q[2], q[3]);
    auto rm = quat.toRotationMatrix();

    const auto m_mat = Eigen::Matrix<double, 6, 6>({
        {m[0][0], m[0][1], m[0][2], m[0][3], m[0][4], m[0][5]},
        {m[1][0], m[1][1], m[1][2], m[1][3], m[1][4], m[1][5]},
        {m[2][0], m[2][1], m[2][2], m[2][3], m[2][4], m[2][5]},
        {m[3][0], m[3][1], m[3][2], m[3][3], m[3][4], m[3][5]},
        {m[4][0], m[4][1], m[4][2], m[4][3], m[4][4], m[4][5]},
        {m[5][0], m[5][1], m[5][2], m[5][3], m[5][4], m[5][5]},
    });

    auto mt = Eigen::Matrix<double, 6, 6>();
    mt.block<3, 6>(0, 0) = rm * m_mat.block<3, 6>(0, 0);
    mt.block<3, 6>(3, 0) = rm * m_mat.block<3, 6>(3, 0);

    rm.transposeInPlace();
    auto mo = Eigen::Matrix<double, 6, 6>();
    mo.block<6, 3>(0, 0) = mt.block<6, 3>(0, 0) * rm;
    mo.block<6, 3>(0, 3) = mt.block<6, 3>(0, 3) * rm;

    return std::array{
        std::array{mo(0, 0), mo(0, 1), mo(0, 2), mo(0, 3), mo(0, 4), mo(0, 5)},
        std::array{mo(1, 0), mo(1, 1), mo(1, 2), mo(1, 3), mo(1, 4), mo(1, 5)},
        std::array{mo(2, 0), mo(2, 1), mo(2, 2), mo(2, 3), mo(2, 4), mo(2, 5)},
        std::array{mo(3, 0), mo(3, 1), mo(3, 2), mo(3, 3), mo(3, 4), mo(3, 5)},
        std::array{mo(4, 0), mo(4, 1), mo(4, 2), mo(4, 3), mo(4, 4), mo(4, 5)},
        std::array{mo(5, 0), mo(5, 1), mo(5, 2), mo(5, 3), mo(5, 4), mo(5, 5)},
    };
}

}  // namespace kynema::math
