#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::state {

/**
 * @brief A Kernel to update the absolute position of each node based on the solver's current
 * state and the initial absolute position.
 */
template <typename DeviceType>
struct UpdateGlobalPosition {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<double* [7]> q;
    ConstView<double* [7]> x0;
    View<double* [7]> x;

    KOKKOS_FUNCTION
    void operator()(int node) const {
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        // Calculate global position
        x(node, 0) = x0(node, 0) + q(node, 0);
        x(node, 1) = x0(node, 1) + q(node, 1);
        x(node, 2) = x0(node, 2) + q(node, 2);

        // Calculate global orientation
        auto node_x_data = Array<double, 4>{};
        const auto node_x = View<double[4]>{node_x_data.data()};
        math::QuaternionCompose(
            subview(q, node, make_pair(3, 7)), subview(x0, node, make_pair(3, 7)), node_x
        );
        x(node, 3) = node_x(0);
        x(node, 4) = node_x(1);
        x(node, 5) = node_x(2);
        x(node, 6) = node_x(3);
    }
};

}  // namespace kynema_fmb::state
