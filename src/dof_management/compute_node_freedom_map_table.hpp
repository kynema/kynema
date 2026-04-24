#pragma once

#include <Kokkos_Core.hpp>

#include "freedom_signature.hpp"
#include "state/state.hpp"

namespace kynema_fmb::dof {

/**
 * @brief A Scanning Kernel which to convert the number of active degrees of freedom per node to
 * a pointer map to the start of their degrees of freedom in a serialized global vector
 */
template <typename DeviceType>
struct ComputeNodeFreedomMapTable {
    typename Kokkos::View<FreedomSignature*, DeviceType>::const_type node_freedom_allocation_table;
    Kokkos::View<size_t*, DeviceType> node_freedom_map_table;

    KOKKOS_FUNCTION
    void operator()(size_t i, size_t& update, bool is_final) const {
        const auto num_dof = count_active_dofs(node_freedom_allocation_table(i));
        update += num_dof;
        if (is_final) {
            node_freedom_map_table(i + 1) = update;
        }
    }
};

/**
 * @brief Compute the node freedom map table, a pointer map to the start of the degrees of freedom
 * ofa given node in a serialized global vector
 *
 * @tparam DeviceType The Kokkos Device where the State object resides
 *
 * @param state A state object with a completed node freedom allocation table
 */
template <typename DeviceType>
inline void compute_node_freedom_map_table(State<DeviceType>& state) {
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    Kokkos::deep_copy(state.node_freedom_map_table, 0UL);
    auto result = 0UL;
    auto scan_range = RangePolicy(0, state.num_system_nodes - 1UL);
    Kokkos::parallel_scan(
        "Compute Node Freedom Map Table", scan_range,
        ComputeNodeFreedomMapTable<DeviceType>{
            state.node_freedom_allocation_table, state.node_freedom_map_table
        },
        result
    );
}

}  // namespace kynema_fmb::dof
