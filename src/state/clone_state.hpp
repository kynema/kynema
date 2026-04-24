#pragma once

#include "state.hpp"

namespace kynema_fmb {

/**
 * @brief Creates a new state object and performs a deep copy of the data in the old one.
 * This is primarily for creating an identical state for snapshotting and rollback in the
 * event that a time step should be performed again.
 *
 * @tparam DeviceType The Kokkos Device where the old and new states live
 * @param old The State to be cloned
 * @return A new State with contents identical to the input
 */
template <typename DeviceType>
inline State<DeviceType> CloneState(const State<DeviceType>& old) {
    using Kokkos::deep_copy;

    auto clone = State<DeviceType>(old.num_system_nodes);
    clone.time_step = old.time_step;
    deep_copy(clone.ID, old.ID);
    deep_copy(clone.node_freedom_allocation_table, old.node_freedom_allocation_table);
    deep_copy(clone.node_freedom_map_table, old.node_freedom_map_table);
    deep_copy(clone.x0, old.x0);
    deep_copy(clone.x, old.x);
    deep_copy(clone.q_delta, old.q_delta);
    deep_copy(clone.q_prev, old.q_prev);
    deep_copy(clone.q, old.q);
    deep_copy(clone.v, old.v);
    deep_copy(clone.vd, old.vd);
    deep_copy(clone.a, old.a);
    deep_copy(clone.f, old.f);
    deep_copy(clone.tangent, old.tangent);
    return clone;
}

}  // namespace kynema_fmb
