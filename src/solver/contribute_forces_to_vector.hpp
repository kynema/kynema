#pragma once

#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which sums the nodal forces into the global RHS vector
 */
template <typename DeviceType>
struct ContributeForcesToVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    ConstView<dof::FreedomSignature*> node_freedom_allocation_table;
    ConstView<size_t*> node_freedom_map_table;
    ConstView<double**> node_loads;
    LeftView<double* [1]> vector;

    KOKKOS_FUNCTION
    void operator()(size_t node) const {
        const auto this_node_freedom_signature = node_freedom_allocation_table(node);
        const auto dof_index = node_freedom_map_table(node);
        if (this_node_freedom_signature == dof::FreedomSignature::NoComponents) {
            return;
        } else if (this_node_freedom_signature == dof::FreedomSignature::AllComponents) {
            vector(dof_index + 0, 0) -= node_loads(node, 0);
            vector(dof_index + 1, 0) -= node_loads(node, 1);
            vector(dof_index + 2, 0) -= node_loads(node, 2);
            vector(dof_index + 3, 0) -= node_loads(node, 3);
            vector(dof_index + 4, 0) -= node_loads(node, 4);
            vector(dof_index + 5, 0) -= node_loads(node, 5);
        } else if (this_node_freedom_signature == dof::FreedomSignature::JustPosition) {
            vector(dof_index + 0, 0) -= node_loads(node, 0);
            vector(dof_index + 1, 0) -= node_loads(node, 1);
            vector(dof_index + 2, 0) -= node_loads(node, 2);
        } else if (this_node_freedom_signature == dof::FreedomSignature::JustRotation) {
            vector(dof_index + 0, 0) -= node_loads(node, 3);
            vector(dof_index + 1, 0) -= node_loads(node, 4);
            vector(dof_index + 2, 0) -= node_loads(node, 5);
        }
    }
};

}  // namespace kynema_fmb::solver
