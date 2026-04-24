#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Kernel for calculating the force applied to the system residual as the result of
 * a revolute joint constraint
 */
template <typename DeviceType>
struct CalculateRevoluteJointForce {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3][3]>& axes, const ConstView<double[7]>& inputs,
        const ConstView<double[7]>& node_u, const View<double[6]>& system_residual_terms
    ) {
        using Kokkos::Array;

        // Initial difference between nodes
        const auto X0_data = Array<double, 3>{axes(0, 0), axes(0, 1), axes(0, 2)};
        const auto R2_data = Array<double, 4>{node_u(3), node_u(4), node_u(5), node_u(6)};
        auto R2_X0_data = Array<double, 3>{};

        const auto X0 = ConstView<double[3]>{X0_data.data()};
        const auto R2 = ConstView<double[4]>{R2_data.data()};
        const auto R2_X0 = View<double[3]>{R2_X0_data.data()};

        //----------------------------------------------------------------------
        // Residual Vector
        //----------------------------------------------------------------------

        // Extract residual rows relevant to this constraint
        math::RotateVectorByQuaternion(R2, X0, R2_X0);

        // Take axis_x and rotate it to right orientation
        system_residual_terms(3) = R2_X0(0) * inputs(0);
        system_residual_terms(4) = R2_X0(1) * inputs(0);
        system_residual_terms(5) = R2_X0(2) * inputs(0);
    }
};

}  // namespace kynema_fmb::constraints
