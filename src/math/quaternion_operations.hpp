#pragma once

#include <array>
#include <cmath>
#include <numbers>

#include <Eigen/Geometry>
#include <Kokkos_Core.hpp>

namespace kynema::math {

/**
 * @brief Converts a 4x1 quaternion to a 3x3 rotation matrix and returns the result
 */
template <typename Quaternion, typename RotationMatrix>
KOKKOS_INLINE_FUNCTION void QuaternionToRotationMatrix(
    const Quaternion& q, const RotationMatrix& R
) {
    R(0, 0) = q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3);
    R(0, 1) = 2. * (q(1) * q(2) - q(0) * q(3));
    R(0, 2) = 2. * (q(1) * q(3) + q(0) * q(2));
    R(1, 0) = 2. * (q(1) * q(2) + q(0) * q(3));
    R(1, 1) = q(0) * q(0) - q(1) * q(1) + q(2) * q(2) - q(3) * q(3);
    R(1, 2) = 2. * (q(2) * q(3) - q(0) * q(1));
    R(2, 0) = 2. * (q(1) * q(3) - q(0) * q(2));
    R(2, 1) = 2. * (q(2) * q(3) + q(0) * q(1));
    R(2, 2) = q(0) * q(0) - q(1) * q(1) - q(2) * q(2) + q(3) * q(3);
}

/**
 * @brief Rotates provided vector by provided *unit* quaternion and returns the result
 */
template <typename Quaternion, typename View1, typename View2>
KOKKOS_INLINE_FUNCTION void RotateVectorByQuaternion(
    const Quaternion& q, const View1& v, const View2& v_rot
) {
    v_rot(0) = (q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3)) * v(0) +
               2. * (q(1) * q(2) - q(0) * q(3)) * v(1) + 2. * (q(1) * q(3) + q(0) * q(2)) * v(2);
    v_rot(1) = 2. * (q(1) * q(2) + q(0) * q(3)) * v(0) +
               (q(0) * q(0) - q(1) * q(1) + q(2) * q(2) - q(3) * q(3)) * v(1) +
               2. * (q(2) * q(3) - q(0) * q(1)) * v(2);
    v_rot(2) = 2. * (q(1) * q(3) - q(0) * q(2)) * v(0) + 2. * (q(2) * q(3) + q(0) * q(1)) * v(1) +
               (q(0) * q(0) - q(1) * q(1) - q(2) * q(2) + q(3) * q(3)) * v(2);
}

/**
 * @brief Computes the derivative of a quaternion and stores the result in a 3x4 matrix
 */
template <typename Quaternion, typename Matrix>
KOKKOS_INLINE_FUNCTION void QuaternionDerivative(const Quaternion& q, const Matrix& m) {
    m(0, 0) = -q(1);
    m(0, 1) = q(0);
    m(0, 2) = -q(3);
    m(0, 3) = q(2);
    m(1, 0) = -q(2);
    m(1, 1) = q(3);
    m(1, 2) = q(0);
    m(1, 3) = -q(1);
    m(2, 0) = -q(3);
    m(2, 1) = -q(2);
    m(2, 2) = q(1);
    m(2, 3) = q(0);
}

/**
 * @brief Computes the inverse of a quaternion
 */
template <typename QuaternionInput, typename QuaternionOutput>
KOKKOS_INLINE_FUNCTION void QuaternionInverse(
    const QuaternionInput& q_in, const QuaternionOutput& q_out
) {
    auto length =
        Kokkos::sqrt(q_in(0) * q_in(0) + q_in(1) * q_in(1) + q_in(2) * q_in(2) + q_in(3) * q_in(3));

    // Inverse of a quaternion is the conjugate divided by the length
    q_out(0) = q_in(0) / length;
    for (auto i = 1; i < 4; ++i) {
        q_out(i) = -q_in(i) / length;
    }
}

/**
 * @brief Composes (i.e. multiplies) two quaternions and stores the result in a third quaternion
 */
template <typename Quaternion1, typename Quaternion2, typename QuaternionN>
KOKKOS_INLINE_FUNCTION void QuaternionCompose(
    const Quaternion1& q1, const Quaternion2& q2, QuaternionN& qn
) {
    qn(0) = q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3);
    qn(1) = q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2);
    qn(2) = q1(0) * q2(2) - q1(1) * q2(3) + q1(2) * q2(0) + q1(3) * q2(1);
    qn(3) = q1(0) * q2(3) + q1(1) * q2(2) - q1(2) * q2(1) + q1(3) * q2(0);
}

/**
 * @brief Returns a 4-D quaternion from provided 3-D rotation vector, i.e. the exponential map
 */
template <typename Vector, typename Quaternion>
KOKKOS_INLINE_FUNCTION void RotationVectorToQuaternion(
    const Vector& phi, const Quaternion& quaternion
) {
    const auto angle = Kokkos::sqrt(phi(0) * phi(0) + phi(1) * phi(1) + phi(2) * phi(2));
    const auto cos_angle = Kokkos::cos(angle / 2.);
    const auto factor = (Kokkos::abs(angle) < 1.e-12) ? 0. : Kokkos::sin(angle / 2.) / angle;

    quaternion(0) = cos_angle;
    for (auto i = 1; i < 4; ++i) {
        quaternion(i) = phi(i - 1) * factor;
    }
}

/**
 * @brief Returns a 3-D rotation vector from provided 4-D quaternion, i.e. the logarithmic map
 */
template <typename Quaternion, typename Vector>
KOKKOS_INLINE_FUNCTION void QuaternionToRotationVector(
    const Quaternion& quaternion, const Vector& phi
) {
    auto theta = 2. * Kokkos::acos(quaternion(0));
    const auto sin_half_theta = std::sqrt(1. - quaternion(0) * quaternion(0));
    if (sin_half_theta > 1e-12) {
        phi(0) = theta * quaternion(1) / sin_half_theta;
        phi(1) = theta * quaternion(2) / sin_half_theta;
        phi(2) = theta * quaternion(3) / sin_half_theta;
    } else {
        phi(0) = 0.;
        phi(1) = 0.;
        phi(2) = 0.;
    }
}

/**
 * @brief Normalizes a quaternion to ensure it is a unit quaternion
 *
 * @details If the length of the quaternion is zero, it returns a default unit quaternion.
 * Otherwise, it normalizes the quaternion and returns the result.
 *
 * @param q The input quaternion as a Kokkos::Array<double, 4>
 * @return Kokkos::Array<double, 4> The normalized quaternion
 */
KOKKOS_INLINE_FUNCTION
Kokkos::Array<double, 4> NormalizeQuaternion(const Kokkos::Array<double, 4>& q) {
    const auto length_squared = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];

    // If the length is 1, our work is done
    if (std::abs(length_squared - 1.) < 1.e-16) {
        return q;
    }

    // If the length of the quaternion is zero, return a default unit quaternion
    if (std::abs(length_squared) < 1.e-16) {
        return Kokkos::Array<double, 4>{1., 0., 0., 0.};
    }

    // Normalize the quaternion
    const auto length = std::sqrt(length_squared);
    auto normalized_quaternion = Kokkos::Array<double, 4>{};
    for (auto k = 0U; k < 4; ++k) {
        normalized_quaternion[k] = q[k] / length;
    }
    return normalized_quaternion;
}

/**
 * @brief Returns a 4-D quaternion from provided tangent vector and twist (degrees) about tangent
 */
inline std::array<double, 4> TangentTwistToQuaternion(
    const std::array<double, 3>& tangent, const double twist
) {
    const auto tan_vec = Eigen::Matrix<double, 3, 1>(tangent.data());
    const auto e1 = tan_vec.normalized();

    const auto temp = Eigen::Matrix<double, 3, 1>::Unit(1);
    const auto plane_axis = (std::abs(e1.dot(temp)) > .9) ? Eigen::Matrix<double, 3, 1>::Unit(0)
                                                          : Eigen::Matrix<double, 3, 1>::Unit(1);
    const auto a = e1.dot(plane_axis);

    const auto e2 = (plane_axis - a * e1).normalized();
    const auto e3 = e1.cross(e2);

    const auto q_tan = Eigen::Quaternion<double>(Eigen::Matrix<double, 3, 3>(
        {{e1(0), e2(0), e3(0)}, {e1(1), e2(1), e3(1)}, {e1(2), e2(2), e3(2)}}
    ));
    const auto twist_rad = twist * std::numbers::pi / 180.;
    const auto q_twist = Eigen::Quaternion<double>(Eigen::AngleAxis<double>(twist_rad, e1));
    const auto q_final = q_twist * q_tan;
    return std::array{q_final.w(), q_final.x(), q_final.y(), q_final.z()};
}

/**
 * @brief Checks if a quaternion is approximately the identity quaternion [1, 0, 0, 0]
 *
 * @param q The quaternion to check
 * @param tolerance The tolerance for the comparison (default: 1e-12)
 * @return true if the quaternion is approximately the identity quaternion, false otherwise
 */
inline bool IsIdentityQuaternion(const std::array<double, 4>& q, double tolerance = 1e-12) {
    return std::abs(q[0] - 1.) <= tolerance && std::abs(q[1]) <= tolerance &&
           std::abs(q[2]) <= tolerance && std::abs(q[3]) <= tolerance;
}

}  // namespace kynema::math
