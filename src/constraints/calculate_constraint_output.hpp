#pragma once

#include <Kokkos_Core.hpp>

#include "calculate_revolute_joint_output.hpp"
#include "constraint_type.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Kernel that calculates the output for a constraints, for use as feedback
 * to controllers
 */
template <typename DeviceType>
struct CalculateConstraintOutput {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<ConstraintType*> type;
    ConstView<size_t*> target_node_index;
    ConstView<double* [3][3]> axes;
    ConstView<double* [7]> node_x0;     // Initial position
    ConstView<double* [7]> node_u;      // Displacement
    ConstView<double* [6]> node_udot;   // Velocity
    ConstView<double* [6]> node_uddot;  // Acceleration
    View<double* [3]> outputs;

    KOKKOS_FUNCTION
    void operator()(int constraint) const {
        if (type(constraint) == ConstraintType::RevoluteJoint) {
            CalculateRevoluteJointOutput<DeviceType>{constraint, target_node_index,
                                                     axes,       node_x0,
                                                     node_u,     node_udot,
                                                     node_uddot, outputs}();
            return;
        }
    }
};

}  // namespace kynema_fmb::constraints
