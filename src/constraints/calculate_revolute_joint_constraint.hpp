#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Kernel for calculating the residual and system gradient for a revolute joint constraint
 */
template <typename DeviceType>
struct CalculateRevoluteJointConstraint {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& X0, const ConstView<double[3][3]>& axes,
        const ConstView<double[7]>& base_node_u, const ConstView<double[7]>& target_node_u,
        const View<double[6]>& residual_terms, const View<double[6][6]>& base_gradient_terms,
        const View<double[6][6]>& target_gradient_terms
    ) {
        using Kokkos::Array;

        const auto u1_data = Array<double, 3>{base_node_u(0), base_node_u(1), base_node_u(2)};
        const auto R1_data =
            Array<double, 4>{base_node_u(3), base_node_u(4), base_node_u(5), base_node_u(6)};
        const auto u2_data = Array<double, 3>{target_node_u(0), target_node_u(1), target_node_u(2)};
        const auto R2_data =
            Array<double, 4>{target_node_u(3), target_node_u(4), target_node_u(5), target_node_u(6)};
        const auto x0_data = Array<double, 3>{axes(0, 0), axes(0, 1), axes(0, 2)};
        const auto y0_data = Array<double, 3>{axes(1, 0), axes(1, 1), axes(1, 2)};
        const auto z0_data = Array<double, 3>{axes(2, 0), axes(2, 1), axes(2, 2)};
        auto R1t_data = Array<double, 4>{};
        auto R1_X0_data = Array<double, 4>{};
        auto x_data = Array<double, 3>{};
        auto y_data = Array<double, 3>{};
        auto z_data = Array<double, 3>{};
        auto xcy_data = Array<double, 3>{};
        auto xcz_data = Array<double, 3>{};

        const auto u1 = ConstView<double[3]>{u1_data.data()};
        const auto R1 = ConstView<double[4]>{R1_data.data()};
        const auto u2 = ConstView<double[3]>{u2_data.data()};
        const auto R2 = ConstView<double[4]>{R2_data.data()};
        const auto x0 = ConstView<double[3]>{x0_data.data()};
        const auto y0 = ConstView<double[3]>{y0_data.data()};
        const auto z0 = ConstView<double[3]>{z0_data.data()};
        const auto R1t = View<double[4]>{R1t_data.data()};
        const auto R1_X0 = View<double[4]>{R1_X0_data.data()};
        const auto x = View<double[3]>{x_data.data()};
        const auto y = View<double[3]>{y_data.data()};
        const auto z = View<double[3]>{z_data.data()};
        const auto xcy = View<double[3]>{xcy_data.data()};
        const auto xcz = View<double[3]>{xcz_data.data()};

        //----------------------------------------------------------------------
        // Residual Vector, Phi
        //----------------------------------------------------------------------

        // Phi(0:3) = u2 + X0 - u1 - R1*X0
        math::QuaternionInverse(R1, R1t);
        math::RotateVectorByQuaternion(R1, X0, R1_X0);
        for (auto component = 0; component < 3; ++component) {
            residual_terms(component) =
                u2(component) + X0(component) - u1(component) - R1_X0(component);
        }

        // Angular residual
        math::RotateVectorByQuaternion(R1, x0, x);
        math::RotateVectorByQuaternion(R2, y0, y);
        math::RotateVectorByQuaternion(R2, z0, z);
        // Phi(3) = dot(R2 * z0_hat, R1 * x0_hat)
        residual_terms(3) = math::DotProduct(z, x);
        // Phi(4) = dot(R2 * y0_hat, R1 * x0_hat)
        residual_terms(4) = math::DotProduct(y, x);

        //----------------------------------------------------------------------
        // Constraint Gradient Matrix, B
        //----------------------------------------------------------------------

        //---------------------------------
        // Target Node
        //---------------------------------
        for (auto component = 0; component < 3; ++component) {
            target_gradient_terms(component, component) = 1.;
        }

        // B(3, 3:6) = -cross(R1 * x0_hat, transpose(R2 * z0_hat))
        math::CrossProduct(x, z, xcz);
        // B(4, 3:6) = -cross(R1 * x0_hat, transpose(R2 * y0_hat))
        math::CrossProduct(x, y, xcy);
        for (auto component = 0; component < 3; ++component) {
            target_gradient_terms(3, component + 3) = -xcz(component);
            target_gradient_terms(4, component + 3) = -xcy(component);
        }
        //---------------------------------
        // Base Node
        //---------------------------------
        // B(0:3,0:3) = -I
        for (auto component = 0; component < 3; ++component) {
            base_gradient_terms(component, component) = -1.;
        }

        // B(3,3:6) = cross(R1 * x0_hat, transpose(R2 * z0_hat))
        // B(4,3:6) = cross(R1 * x0_hat, transpose(R2 * y0_hat))
        for (auto component = 0; component < 3; ++component) {
            base_gradient_terms(3, component + 3) = xcz(component);
            base_gradient_terms(4, component + 3) = xcy(component);
        }
    }
};
}  // namespace kynema_fmb::constraints
