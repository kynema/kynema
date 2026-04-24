#pragma once

#include <ranges>
#include <span>
#include <vector>

#include "node.hpp"
#include "state/state.hpp"
#include "state/update_global_position.hpp"

namespace kynema_fmb::model {

template <typename DeviceType>
inline void CopyNodesToState(State<DeviceType>& state, std::span<const Node> nodes) {
    using Kokkos::create_mirror_view;
    using Kokkos::deep_copy;
    using Kokkos::WithoutInitializing;
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;

    auto host_id = create_mirror_view(WithoutInitializing, state.ID);
    auto host_x0 = create_mirror_view(WithoutInitializing, state.x0);
    auto host_q = create_mirror_view(WithoutInitializing, state.q);
    auto host_v = create_mirror_view(WithoutInitializing, state.v);
    auto host_vd = create_mirror_view(WithoutInitializing, state.vd);

    for (auto node_index : std::views::iota(0U, nodes.size())) {
        const auto& node = nodes[node_index];
        for (auto component : std::views::iota(0U, 7U)) {
            host_x0(node_index, component) = node.x0[component];
            host_q(node_index, component) = node.u[component];
        }
        for (auto component : std::views::iota(0U, 6U)) {
            host_v(node_index, component) = node.v[component];
            host_vd(node_index, component) = node.vd[component];
        }
        host_id(node_index) = node.id;
    }

    // Copy data to host
    deep_copy(state.ID, host_id);
    deep_copy(state.x0, host_x0);
    deep_copy(state.q, host_q);
    deep_copy(state.v, host_v);
    deep_copy(state.vd, host_vd);
    deep_copy(state.a, state.vd);  // initialize algorithmic acceleration to acceleration
    deep_copy(state.f, 0.);

    // Set previous state to current state
    deep_copy(state.q_prev, state.q);

    // Calculate current global position from initial position and displacement
    Kokkos::parallel_for(
        "UpdateGlobalPosition", RangePolicy(0UL, state.num_system_nodes),
        state::UpdateGlobalPosition<DeviceType>{
            state.q,
            state.x0,
            state.x,
        }
    );
}

}  // namespace kynema_fmb::model
