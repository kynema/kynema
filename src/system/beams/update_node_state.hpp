#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct UpdateNodeStateElement {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element{};
    ConstView<size_t**> node_state_indices;
    View<double* [7]> node_u;
    View<double* [6]> node_u_dot;
    View<double* [6]> node_u_ddot;

    ConstView<double* [7]> Q;
    ConstView<double* [6]> V;
    ConstView<double* [6]> A;

    KOKKOS_FUNCTION
    void operator()(const size_t node) const {
        const auto index = node_state_indices(element, node);
        for (auto component = 0U; component < 7U; ++component) {
            node_u(node, component) = Q(index, component);
        }
        for (auto component = 0U; component < 6U; ++component) {
            node_u_dot(node, component) = V(index, component);
            node_u_ddot(node, component) = A(index, component);
        }
    }
};

template <typename DeviceType>
struct UpdateNodeState {
    using TeamPolicy = typename Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<double* [7]> Q;
    ConstView<double* [6]> V;
    ConstView<double* [6]> A;
    ConstView<size_t**> node_state_indices;
    ConstView<size_t*> num_nodes_per_element;
    View<double** [7]> node_u_;
    View<double** [6]> node_u_dot_;
    View<double** [6]> node_u_ddot_;

    KOKKOS_FUNCTION
    void operator()(member_type member) const {
        const auto element = static_cast<size_t>(member.league_rank());
        const auto num_nodes = num_nodes_per_element(element);
        const auto node_range = Kokkos::TeamThreadRange(member, num_nodes);

        const auto node_u = View<double* [7]>(member.team_scratch(1), num_nodes);
        const auto node_u_dot = View<double* [6]>(member.team_scratch(1), num_nodes);
        const auto node_u_ddot = View<double* [6]>(member.team_scratch(1), num_nodes);

        const auto node_state_updater = beams::UpdateNodeStateElement<DeviceType>{
            element, node_state_indices, node_u, node_u_dot, node_u_ddot, Q, V, A
        };
        Kokkos::parallel_for(node_range, node_state_updater);

        using CopyMatrix = KokkosBatched::TeamVectorCopy<member_type>;

        CopyMatrix::invoke(
            member, node_u, Kokkos::subview(node_u_, element, Kokkos::ALL, Kokkos::ALL)
        );
        CopyMatrix::invoke(
            member, node_u_dot, Kokkos::subview(node_u_dot_, element, Kokkos::ALL, Kokkos::ALL)
        );
        CopyMatrix::invoke(
            member, node_u_ddot, Kokkos::subview(node_u_ddot_, element, Kokkos::ALL, Kokkos::ALL)
        );
    }
};

}  // namespace kynema_fmb::beams
