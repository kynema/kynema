#pragma once

#include <Kokkos_Core.hpp>

#include "elements/beams/beams.hpp"
#include "elements/masses/masses.hpp"
#include "elements/springs/springs.hpp"

namespace kynema_fmb {

/**
 * @brief A container providing handle to all structural elements present in the model
 *
 * @details This class serves as a container for wrapping all structural elements
 *          in the model, including beams, point masses, and springs, effectively providing
 *          the entire mesh as a single object to be used by the gen-alpha solver
 *          logic.
 */
template <typename DeviceType>
struct Elements {
    Beams<DeviceType> beams;
    Masses<DeviceType> masses;
    Springs<DeviceType> springs;

    Elements() : beams(0U, 0U, 0U), masses(0U), springs(0U) {}

    Elements(Beams<DeviceType> b, Masses<DeviceType> m, Springs<DeviceType> s)
        : beams(std::move(b)), masses(std::move(m)), springs(std::move(s)) {}

    /// Returns the total number of elements across all types in the system
    [[nodiscard]] size_t NumElementsInSystem() const {
        return beams.num_elems + masses.num_elems + springs.num_elems;
    }

    /**
     * @brief Returns the number of nodes per element for each element in the system
     *
     * @return Kokkos::View<size_t*> A 1D view containing the number of nodes for each element,
     *         concatenated in the order: [beams] -> [masses] -> [springs]
     */
    [[nodiscard]] Kokkos::View<size_t*, DeviceType> NumberOfNodesPerElement() const {
        using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
        const auto result = Kokkos::View<size_t*, DeviceType>(
            Kokkos::view_alloc("num_nodes_per_element", Kokkos::WithoutInitializing),
            NumElementsInSystem()
        );

        // Beams
        auto beams_num_nodes_per_element = beams.num_nodes_per_element;
        auto beams_range = RangePolicy(0, beams.num_elems);
        Kokkos::parallel_for(
            beams_range,
            KOKKOS_LAMBDA(size_t element) { result(element) = beams_num_nodes_per_element(element); }
        );

        // Masses
        auto beams_offset = beams.num_elems;
        auto masses_num_nodes_per_element = masses.num_nodes_per_element;
        auto masses_range = RangePolicy(0, masses.num_elems);
        Kokkos::parallel_for(
            masses_range, KOKKOS_LAMBDA(size_t element) {
                result(element + beams_offset) = masses_num_nodes_per_element(element);
            }
        );

        // Springs
        auto beams_and_masses_offset = beams_offset + masses.num_elems;
        auto springs_num_nodes_per_element = springs.num_nodes_per_element;
        auto springs_range = RangePolicy(0, springs.num_elems);
        Kokkos::parallel_for(
            springs_range, KOKKOS_LAMBDA(size_t element) {
                result(element + beams_and_masses_offset) = springs_num_nodes_per_element(element);
            }
        );

        return result;
    }

    /**
     * @brief Returns the state indices for each node of each element in the system
     *
     * @return Kokkos::View<size_t**> A 2D view containing the state indices for each node,
     *         concatenated in the order: [beams] -> [masses] -> [springs]
     */
    [[nodiscard]] Kokkos::View<size_t**, DeviceType> NodeStateIndices() const {
        using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
        const auto max_nodes = std::max(beams.max_elem_nodes, springs.num_elems > 0 ? 2UL : 1UL);
        const auto result = Kokkos::View<size_t**, DeviceType>(
            Kokkos::view_alloc("node_state_indices", Kokkos::WithoutInitializing),
            NumElementsInSystem(), max_nodes
        );

        // Beams
        auto beams_num_nodes_per_element = beams.num_nodes_per_element;
        auto beams_node_state_indices = beams.node_state_indices;
        auto beams_range = RangePolicy(0, beams.num_elems);
        Kokkos::parallel_for(
            beams_range, KOKKOS_LAMBDA(size_t element) {
                const auto num_nodes = beams_num_nodes_per_element(element);
                for (auto node = 0U; node < num_nodes; ++node) {
                    result(element, node) = beams_node_state_indices(element, node);
                }
            }
        );

        // Masses
        const auto beams_offset = beams.num_elems;
        auto masses_state_indices = masses.state_indices;
        auto masses_range = RangePolicy(0, masses.num_elems);
        Kokkos::parallel_for(
            masses_range, KOKKOS_LAMBDA(size_t element) {
                // Masses always have one node per element
                result(element + beams_offset, 0) = masses_state_indices(element);
            }
        );

        // Springs
        auto beams_and_masses_offset = beams_offset + masses.num_elems;
        auto springs_node_state_indices = springs.node_state_indices;
        auto springs_range = RangePolicy(0, springs.num_elems);
        Kokkos::parallel_for(
            springs_range, KOKKOS_LAMBDA(size_t element) {
                // Springs always have two nodes per element
                result(element + beams_and_masses_offset, 0) =
                    springs_node_state_indices(element, 0);
                result(element + beams_and_masses_offset, 1) =
                    springs_node_state_indices(element, 1);
            }
        );

        return result;
    }
};

}  // namespace kynema_fmb
