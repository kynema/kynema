#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Kernel for calculating the output for a revolute joint constraint for feedback to
 * controllers
 */
template <typename DeviceType>
struct CalculateRevoluteJointOutput {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    int constraint;
    ConstView<size_t*> target_node_index;
    ConstView<double* [3][3]> axes;
    ConstView<double* [7]> node_x0;     // Initial position
    ConstView<double* [7]> node_u;      // Displacement
    ConstView<double* [6]> node_udot;   // Velocity
    ConstView<double* [6]> node_uddot;  // Acceleration
    View<double* [3]> outputs;

    KOKKOS_FUNCTION
    void operator()() const {
        using Kokkos::Array;

        // Axis of rotation unit vector
        const auto joint_axis0_data =
            Array<double, 3>{axes(constraint, 0, 0), axes(constraint, 0, 1), axes(constraint, 0, 2)};
        const auto joint_axis0 = ConstView<double[3]>{joint_axis0_data.data()};

        // Target node index
        auto node = target_node_index(constraint);

        // Target node initial rotation
        const auto R0_data =
            Array<double, 4>{node_u(node, 3), node_u(node, 4), node_u(node, 5), node_u(node, 6)};
        const auto R0 = ConstView<double[4]>{R0_data.data()};

        // Target node rotational displacement
        const auto R_data =
            Array<double, 4>{node_u(node, 3), node_u(node, 4), node_u(node, 5), node_u(node, 6)};
        const auto R = ConstView<double[4]>{R_data.data()};

        // Calculate current orientation
        auto RR0_data = Array<double, 4>{};
        auto RR0 = View<double[4]>(RR0_data.data());
        math::QuaternionCompose(R, R0, RR0);

        // Calculate rotational displacement as rotation vector
        auto RotVec_data = Array<double, 3>{};
        auto RotVec = View<double[3]>(RotVec_data.data());
        math::QuaternionToRotationVector(R, RotVec);

        // Target node rotational velocity vector
        const auto omega_data =
            Array<double, 3>{node_udot(node, 3), node_udot(node, 4), node_udot(node, 5)};
        auto omega = ConstView<double[3]>{omega_data.data()};

        // Target node rotational acceleration vector
        const auto omega_dot_data =
            Array<double, 3>{node_uddot(node, 3), node_uddot(node, 4), node_uddot(node, 5)};
        auto omega_dot = ConstView<double[3]>{omega_dot_data.data()};

        // Calculate joint axis in current configuration
        auto joint_axis_data = Array<double, 3>{};
        auto joint_axis = View<double[3]>{joint_axis_data.data()};
        math::RotateVectorByQuaternion(R, joint_axis0, joint_axis);

        // Calculate rotation about shaft axis
        auto angular_rotation = math::DotProduct(joint_axis, RotVec);

        // Angular velocity about joint axis (rad/s)
        auto angular_velocity = math::DotProduct(joint_axis, omega);

        // Angular acceleration about joint axis (rad/s)
        auto angular_acceleration = math::DotProduct(joint_axis, omega_dot);

        // Save outputs
        outputs(constraint, 0) = angular_rotation;
        outputs(constraint, 1) = angular_velocity;
        outputs(constraint, 2) = angular_acceleration;
    }
};

}  // namespace kynema_fmb::constraints
