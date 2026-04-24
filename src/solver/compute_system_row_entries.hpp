#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief Kernel to compute the elements' contribution to the row pointers of the CRS matrix
 */
template <typename RowPtrType>
struct ComputeSystemRowEntries {
    using ValueType = typename RowPtrType::value_type;
    using DeviceType = typename RowPtrType::device_type;
    template <typename value_type>
    using View = Kokkos::View<value_type, DeviceType>;
    template <typename value_type>
    using ConstView = typename View<value_type>::const_type;

    ConstView<size_t*> active_dofs;
    ConstView<size_t*> node_freedom_map_table;
    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t**> node_state_indices;
    ConstView<size_t*> base_active_dofs;
    ConstView<size_t*> target_active_dofs;
    ConstView<size_t* [6]> base_node_freedom_table;
    ConstView<size_t* [6]> target_node_freedom_table;
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    RowPtrType row_entries;

    KOKKOS_FUNCTION
    bool ElementContainsNode(size_t element, size_t node) const {
        const auto num_nodes = num_nodes_per_element(element);
        for (auto node_index = 0U; node_index < num_nodes; ++node_index) {
            if (node_state_indices(element, node_index) == node) {
                return true;
            }
        }
        return false;
    }

    KOKKOS_FUNCTION
    bool BaseContainsNode(size_t constraint, size_t dof_index) const {
        return dof_index == base_node_freedom_table(constraint, 0) &&
               base_active_dofs(constraint) != 0UL;
    }

    KOKKOS_FUNCTION
    bool TargetContainsNode(size_t constraint, size_t dof_index) const {
        return dof_index == target_node_freedom_table(constraint, 0) &&
               target_active_dofs(constraint) != 0UL;
    }

    KOKKOS_FUNCTION
    bool ConstraintContainsNode(size_t constraint, size_t dof_index) const {
        return BaseContainsNode(constraint, dof_index) || TargetContainsNode(constraint, dof_index);
    }

    KOKKOS_FUNCTION
    ValueType ComputeEntriesInConstraint(size_t constraint) const {
        return static_cast<ValueType>(row_range(constraint).second - row_range(constraint).first);
    }

    KOKKOS_FUNCTION
    ValueType ComputeEntriesInElement(size_t element) const {
        const auto num_nodes = num_nodes_per_element(element);
        auto entries = ValueType{};
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto node_state_index = node_state_indices(element, node);
            entries += static_cast<ValueType>(active_dofs(node_state_index));
        }
        return entries;
    }

    KOKKOS_FUNCTION
    void operator()(size_t node) const {
        const auto num_dof = active_dofs(node);
        const auto num_elements = num_nodes_per_element.extent(0);
        const auto num_constraints = row_range.extent(0);
        const auto dof_index = node_freedom_map_table(node);

        auto entries_in_row = static_cast<ValueType>(num_dof);

        for (auto element = 0U; element < num_elements; ++element) {
            if (!ElementContainsNode(element, node)) {
                continue;
            }
            const auto entries_in_element = ComputeEntriesInElement(element);
            entries_in_row += entries_in_element - static_cast<ValueType>(num_dof);
        }

        for (auto constraint = 0U; constraint < num_constraints; ++constraint) {
            if (!ConstraintContainsNode(constraint, dof_index)) {
                continue;
            }
            const auto entries_in_constraint = ComputeEntriesInConstraint(constraint);
            entries_in_row += entries_in_constraint;
        }

        for (auto component = 0U; component < num_dof; ++component) {
            row_entries(dof_index + component) = entries_in_row;
        }
    }
};

}  // namespace kynema_fmb::solver
