#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which sums the residual contributions computed at each node in a beam
 * into the correct location of the global RHS vector.
 */
template <typename DeviceType>
struct ContributeBeamsToVector {
    using TeamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t** [6]> element_freedom_table;
    ConstView<double** [6]> elements;
    LeftView<double* [1]> vector;

    KOKKOS_FUNCTION
    void operator()(member_type member) const {
        const auto element = static_cast<size_t>(member.league_rank());
        const auto num_nodes = num_nodes_per_element(element);
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;

        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, num_nodes * 6U), [&](size_t node_component) {
                const auto node = node_component % num_nodes;
                const auto component = node_component / num_nodes;
                const auto entry = element_freedom_table(element, node, component);
                if constexpr (force_atomic) {
                    Kokkos::atomic_add(&vector(entry, 0), elements(element, node, component));
                } else {
                    vector(entry, 0) += elements(element, node, component);
                }
            }
        );
    }
};

}  // namespace kynema_fmb::solver
