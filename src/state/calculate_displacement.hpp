#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::state {

/**
 * @brief A Kernel for applying the computed change in state (displacement) to the previous state
 * to get an estimate for the current state.
 */
template <typename DeviceType>
struct CalculateDisplacement {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    double h;
    ConstView<double* [6]> q_delta;
    ConstView<double* [7]> q_prev;
    View<double* [7]> q;

    KOKKOS_FUNCTION
    void operator()(int node) const {
        using Kokkos::Array;

        for (auto component = 0; component < 3; ++component) {
            q(node, component) = q_prev(node, component) + h * q_delta(node, component);
        }

        auto delta_data =
            Array<double, 3>{h * q_delta(node, 3), h * q_delta(node, 4), h * q_delta(node, 5)};
        auto delta = View<double[3]>{delta_data.data()};

        auto quat_delta_data = Array<double, 4>{};
        auto quat_delta = View<double[4]>{quat_delta_data.data()};
        math::RotationVectorToQuaternion(delta, quat_delta);
        auto quat_prev_data =
            Array<double, 4>{q_prev(node, 3), q_prev(node, 4), q_prev(node, 5), q_prev(node, 6)};
        auto quat_prev = View<double[4]>{quat_prev_data.data()};
        auto quat_new_data = Array<double, 4>{};
        auto quat_new = View<double[4]>{quat_new_data.data()};
        math::QuaternionCompose(quat_delta, quat_prev, quat_new);

        for (auto component = 0; component < 4; ++component) {
            q(node, component + 3) = quat_new(component);
        }
    }
};

}  // namespace kynema_fmb::state
