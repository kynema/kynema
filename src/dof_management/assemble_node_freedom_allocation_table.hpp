#pragma once

#include <Kokkos_Core.hpp>

#include "constraints/constraints.hpp"
#include "elements/elements.hpp"
#include "freedom_signature.hpp"
#include "state/state.hpp"

namespace kynema_fmb::dof {

/**
 * @brief A Kernel for applying a Beam element's freedom signature to all nodes it contains
 */
template <typename DeviceType>
struct AssembleNodeFreedomMapTable_Beams {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t**> node_state_indices;
    ConstView<FreedomSignature**> element_freedom_signature;
    View<FreedomSignature*> node_freedom_allocation_table;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        const auto num_nodes = num_nodes_per_element(element);
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto node_index = node_state_indices(element, node);
            Kokkos::atomic_or(
                &node_freedom_allocation_table(node_index), element_freedom_signature(element, node)
            );
        }
    }
};

/**
 * @brief A Kernel for applying a mass element's freedom signature to its node
 */
template <typename DeviceType>
struct AssembleNodeFreedomMapTable_Masses {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> node_state_indices;
    ConstView<FreedomSignature*> element_freedom_signature;
    Kokkos::View<FreedomSignature*, DeviceType> node_freedom_allocation_table;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        // Masses always have one node per element
        const auto node_index = node_state_indices(element);
        Kokkos::atomic_or(
            &node_freedom_allocation_table(node_index), element_freedom_signature(element)
        );
    }
};

/**
 * @brief A Kernel for applying a spring element's freedom signature to both of its nodes
 */
template <typename DeviceType>
struct AssembleNodeFreedomMapTable_Springs {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t* [2]> node_state_indices;
    ConstView<FreedomSignature* [2]> element_freedom_signature;
    View<FreedomSignature*> node_freedom_allocation_table;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        // Springs always have two nodes per element
        for (auto node = 0U; node < 2U; ++node) {
            const auto node_index = node_state_indices(element, node);
            Kokkos::atomic_or(
                &node_freedom_allocation_table(node_index), element_freedom_signature(element, node)
            );
        }
    }
};

/**
 * @brief A Kernel for applying a constraint's freedom signature to its base and target nodes
 */
template <typename DeviceType>
struct AssembleNodeFreedomMapTable_Constraints {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<constraints::ConstraintType*> type;
    ConstView<size_t*> target_node_index;
    ConstView<size_t*> base_node_index;
    ConstView<FreedomSignature*> target_node_freedom_signature;
    ConstView<FreedomSignature*> base_node_freedom_signature;
    View<FreedomSignature*> node_freedom_allocation_table;

    KOKKOS_FUNCTION
    void operator()(size_t node) const {
        {
            const auto node_index = target_node_index(node);
            Kokkos::atomic_or(
                &node_freedom_allocation_table(node_index), target_node_freedom_signature(node)
            );
        }

        if (GetNumberOfNodes(type(node)) == 2U) {
            const auto node_index = base_node_index(node);
            Kokkos::atomic_or(
                &node_freedom_allocation_table(node_index), base_node_freedom_signature(node)
            );
        }
    }
};

/**
 * @brief Creates the node freedom allocation table in state based on the connectivities
 * defined in the elements and constraints structures
 *
 * @details The node freedom allocation table defines which degrees of freedom are defined
 * on which nodes by collecting this information from all elements and constraints which
 * contain this node.  The most common action using this information is to query how many
 * degrees of freedom are defined, so this information is also computed.
 *
 * @tparam DeviceType The Kokkos device defining where Kynema's structures reside
 *
 * @param state A structure containing the node freedom allocation table to be filled
 * @param elements A structure containing all of the elements and their connectivities
 * @param constraints A structure containing all of the constraints and their connectivities
 */
template <typename DeviceType>
inline void assemble_node_freedom_allocation_table(
    State<DeviceType>& state, const Elements<DeviceType>& elements,
    const Constraints<DeviceType>& constraints
) {
    using Kokkos::parallel_for;
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    Kokkos::deep_copy(state.node_freedom_allocation_table, FreedomSignature::NoComponents);

    auto beams_range = RangePolicy(0, elements.beams.num_elems);

    parallel_for(
        "AssembleNodeFreedomMapTable_Beams", beams_range,
        AssembleNodeFreedomMapTable_Beams<DeviceType>{
            elements.beams.num_nodes_per_element, elements.beams.node_state_indices,
            elements.beams.element_freedom_signature, state.node_freedom_allocation_table
        }
    );
    auto masses_range = RangePolicy(0, elements.masses.num_elems);
    parallel_for(
        "AssembleNodeFreedomMapTable_Masses", masses_range,
        AssembleNodeFreedomMapTable_Masses<DeviceType>{
            elements.masses.state_indices, elements.masses.element_freedom_signature,
            state.node_freedom_allocation_table
        }
    );
    auto springs_range = RangePolicy(0, elements.springs.num_elems);
    parallel_for(
        "AssembleNodeFreedomMapTable_Springs", springs_range,
        AssembleNodeFreedomMapTable_Springs<DeviceType>{
            elements.springs.node_state_indices, elements.springs.element_freedom_signature,
            state.node_freedom_allocation_table
        }
    );
    auto constraints_range = RangePolicy(0, constraints.num_constraints);
    parallel_for(
        "AssembleNodeFreedomMapTable_Constraints", constraints_range,
        AssembleNodeFreedomMapTable_Constraints<DeviceType>{
            constraints.type, constraints.target_node_index, constraints.base_node_index,
            constraints.target_node_freedom_signature, constraints.base_node_freedom_signature,
            state.node_freedom_allocation_table
        }
    );

    const auto active_dofs = state.active_dofs;
    const auto node_freedom_allocation_table = state.node_freedom_allocation_table;
    auto system_range = RangePolicy(0, state.num_system_nodes);
    parallel_for(
        "ComputeActiveDofs", system_range,
        KOKKOS_LAMBDA(size_t node) {
            active_dofs(node) = count_active_dofs(node_freedom_allocation_table(node));
        }
    );
}

}  // namespace kynema_fmb::dof
