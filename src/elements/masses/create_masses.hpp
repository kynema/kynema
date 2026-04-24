#pragma once
#include <span>

#include "masses.hpp"
#include "masses_input.hpp"
#include "model/node.hpp"

namespace kynema_fmb {

/**
 * @brief Creates a masses data structure and initializes its data
 *
 * @tparam DeviceData The Kokkos Device where the newly created Masses structure will reside
 *
 * @param masses_input A MassesInput object defining the mass elements
 * @param nodes A vector defining all of the nodes in the problem
 * @return A fully initialized Masses data structure
 *
 */
template <typename DeviceType>
inline Masses<DeviceType> CreateMasses(
    const MassesInput& masses_input, std::span<const Node> nodes
) {
    using Kokkos::create_mirror_view;
    using Kokkos::deep_copy;
    using Kokkos::WithoutInitializing;

    Masses<DeviceType> masses(masses_input.NumElements());

    auto host_gravity = create_mirror_view(WithoutInitializing, masses.gravity);
    auto host_state_indices = create_mirror_view(WithoutInitializing, masses.state_indices);
    auto host_x0 = create_mirror_view(WithoutInitializing, masses.node_x0);
    auto host_Mstar = create_mirror_view(WithoutInitializing, masses.qp_Mstar);

    host_gravity(0) = masses_input.gravity[0];
    host_gravity(1) = masses_input.gravity[1];
    host_gravity(2) = masses_input.gravity[2];

    // Populate element data - x0 and Mstar
    for (auto element : std::views::iota(0U, masses_input.NumElements())) {
        host_state_indices(element) = static_cast<size_t>(masses_input.elements[element].node_id);
        const auto& elem = masses_input.elements[element];
        for (auto component : std::views::iota(0U, 7U)) {
            host_x0(element, component) = nodes[elem.node_id].x0[component];
        }
        for (auto component_1 : std::views::iota(0U, 6U)) {
            for (auto component_2 : std::views::iota(0U, 6U)) {
                host_Mstar(element, component_1, component_2) =
                    elem.M_star[component_1][component_2];
            }
        }
    }

    deep_copy(masses.gravity, host_gravity);
    deep_copy(masses.state_indices, host_state_indices);
    deep_copy(masses.node_x0, host_x0);
    deep_copy(masses.qp_Mstar, host_Mstar);

    return masses;
}

}  // namespace kynema_fmb
