#pragma once

#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::state {

/**
 * @brief A Kernel to update the velocity, acceleration, and change in state at a node for
 * a dynamic problem
 */
template <typename DeviceType>
struct UpdateDynamicPrediction {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    double h;
    double beta_prime;
    double gamma_prime;
    ConstView<dof::FreedomSignature*> node_freedom_allocation_table;
    ConstView<size_t*> node_freedom_map_table;
    ConstLeftView<double* [1]> x_delta;
    View<double* [6]> q_delta;
    View<double* [6]> v;
    View<double* [6]> vd;

    KOKKOS_FUNCTION
    void operator()(size_t node) const {
        const auto num_dof = count_active_dofs(node_freedom_allocation_table(node));
        const auto first_dof = node_freedom_map_table(node);
        for (auto component = 0U; component < num_dof; ++component) {
            const auto delta = x_delta(first_dof + component, 0);
            q_delta(node, component) += delta / h;
            v(node, component) += gamma_prime * delta;
            vd(node, component) += beta_prime * delta;
        }
    }
};

}  // namespace kynema_fmb::state
