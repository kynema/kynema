#pragma once

#include <Kokkos_Core.hpp>

#include "constraints/constraints.hpp"
#include "state/state.hpp"

namespace kynema_fmb::dof {

/**
 * @brief A Kernel that creates the node freedom tables for each the target and base nodes
 * for a given constrain.
 */
template <typename DeviceType>
struct CreateConstraintFreedomTable {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<constraints::ConstraintType*> type;
    ConstView<size_t*> target_node_index;
    ConstView<size_t*> base_node_index;
    ConstView<size_t*> active_dofs;
    ConstView<size_t*> node_freedom_map_table;
    View<size_t* [6]> target_node_freedom_table;
    View<size_t* [6]> base_node_freedom_table;

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        {
            const auto node_index = target_node_index(constraint);
            const auto num_active_dofs = active_dofs(node_index);
            for (auto component = 0U; component < num_active_dofs; ++component) {
                target_node_freedom_table(constraint, component) =
                    node_freedom_map_table(node_index) + component;
            }
        }
        if (GetNumberOfNodes(type(constraint)) == 2U) {
            const auto node_index = base_node_index(constraint);
            const auto num_active_dofs = active_dofs(node_index);
            for (auto component = 0U; component < num_active_dofs; ++component) {
                base_node_freedom_table(constraint, component) =
                    node_freedom_map_table(node_index) + component;
            }
        } else {
            for (auto component = 0U; component < 6U; ++component) {
                base_node_freedom_table(constraint, component) = 0UL;
            }
        }
    }
};

/**
 * @brief Creates node freedom tables for each the target and base nodes in the constraints
 *
 * @details The constraint freedom table maps each degree of freedom for each node of a constraint
 * to its global degree of freedom number
 *
 * @tparam DeviceType The Kokkos Device where constraints and state reside
 *
 * @param constraints The Constraints object used to create state's node freedom map table
 * @param state A State object with a completed node freedom map table
 */
template <typename DeviceType>
inline void create_constraint_freedom_table(
    Constraints<DeviceType>& constraints, const State<DeviceType>& state
) {
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    auto constraints_range = RangePolicy(0, constraints.num_constraints);

    Kokkos::parallel_for(
        "Create Constraint Node Freedom Table", constraints_range,
        CreateConstraintFreedomTable<DeviceType>{
            constraints.type, constraints.target_node_index, constraints.base_node_index,
            state.active_dofs, state.node_freedom_map_table, constraints.target_node_freedom_table,
            constraints.base_node_freedom_table
        }
    );
}

}  // namespace kynema_fmb::dof
