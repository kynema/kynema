#pragma once

#include <Kokkos_Core.hpp>

#include "elements/elements.hpp"
#include "state/state.hpp"

namespace kynema_fmb::dof {

/**
 * @brief A Kernel that creates the element freedom table which maps each degree of freedom
 * on the beam element to its global component number.
 */
template <typename DeviceType>
struct CreateElementFreedomTable_Beams {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t**> node_state_indices;
    ConstView<size_t*> node_freedom_map_table;
    View<size_t***> element_freedom_table;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        for (auto node = 0U; node < num_nodes_per_element(element); ++node) {
            const auto node_index = node_state_indices(element, node);
            for (auto component = 0U; component < 6U; ++component) {
                element_freedom_table(element, node, component) =
                    node_freedom_map_table(node_index) + component;
            }
        }
    }
};

/**
 * @brief A Kernel that creates the element freedom table which maps each degree of freedom
 * on the mass element to its global component number.
 */
template <typename DeviceType>
struct CreateElementFreedomTable_Masses {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> node_state_indices;
    ConstView<size_t*> node_freedom_map_table;
    View<size_t**> element_freedom_table;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        // Masses always have one node per element
        const auto node_index = node_state_indices(element);
        for (auto component = 0U; component < 6U; ++component) {
            element_freedom_table(element, component) =
                node_freedom_map_table(node_index) + component;
        }
    }
};

/**
 * @brief A Kernel that creates the element freedom table which maps each degree of freedom
 * on the spring element to its global component number.
 */
template <typename DeviceType>
struct CreateElementFreedomTable_Springs {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t* [2]> node_state_indices;
    ConstView<size_t*> node_freedom_map_table;
    View<size_t* [2][3]> element_freedom_table;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        // Springs always have two nodes per element
        for (auto node = 0U; node < 2U; ++node) {
            const auto node_index = node_state_indices(element, node);
            // Springs only have translational DOFs
            for (auto component = 0U; component < 3U; ++component) {
                element_freedom_table(element, node, component) =
                    node_freedom_map_table(node_index) + component;
            }
        }
    }
};

/**
 * @brief Creates the element freedom tables for all of the elements in the system
 *
 * @details The element freedom table for each element type maps the degrees of freedom
 * for that element to their global degree of freedom number
 *
 * @tparam DeviceType The Kokkos Device where elements and state reside
 *
 * @param elements the Elements object used to create state's node freedom map table
 * @param state A State object with a completed node freedom map table
 */
template <typename DeviceType>
inline void create_element_freedom_table(
    Elements<DeviceType>& elements, const State<DeviceType>& state
) {
    using Kokkos::parallel_for;
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    auto beams_range = RangePolicy(0, elements.beams.num_elems);

    parallel_for(
        "CreateElementFreedomTable_Beams", beams_range,
        CreateElementFreedomTable_Beams<DeviceType>{
            elements.beams.num_nodes_per_element, elements.beams.node_state_indices,
            state.node_freedom_map_table, elements.beams.element_freedom_table
        }
    );
    auto masses_range = RangePolicy(0, elements.masses.num_elems);
    parallel_for(
        "CreateElementFreedomTable_Masses", masses_range,
        CreateElementFreedomTable_Masses<DeviceType>{
            elements.masses.state_indices, state.node_freedom_map_table,
            elements.masses.element_freedom_table
        }
    );
    auto springs_range = RangePolicy(0, elements.springs.num_elems);
    parallel_for(
        "CreateElementFreedomTable_Springs", springs_range,
        CreateElementFreedomTable_Springs<DeviceType>{
            elements.springs.num_nodes_per_element, elements.springs.node_state_indices,
            state.node_freedom_map_table, elements.springs.element_freedom_table
        }
    );
}

}  // namespace kynema_fmb::dof
