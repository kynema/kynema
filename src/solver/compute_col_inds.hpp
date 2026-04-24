#pragma once

#include <Kokkos_Core.hpp>

#include "compute_constraints_col_inds.hpp"
#include "compute_system_col_inds.hpp"

namespace kynema_fmb::solver {

/**
 * @brief The top level function object for computing the column indicies for the CRS matrix
 * to be solved at each nonlinear iteration
 */
template <typename RowPtrType, typename IndicesType>
struct ComputeColInds {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, typename RowPtrType::device_type>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    [[nodiscard]] static IndicesType invoke(
        typename RowPtrType::value_type num_non_zero, size_t num_system_dofs,
        const ConstView<size_t*>& active_dofs, const ConstView<size_t*>& node_freedom_map_table,
        const ConstView<size_t*>& num_nodes_per_element,
        const ConstView<size_t**>& node_state_indices, const ConstView<size_t*>& base_active_dofs,
        const ConstView<size_t*>& target_active_dofs,
        const ConstView<size_t* [6]>& base_node_freedom_table,
        const ConstView<size_t* [6]>& target_node_freedom_table,
        const ConstView<Kokkos::pair<size_t, size_t>*>& row_range,
        const typename RowPtrType::const_type& row_ptrs
    ) {
        using RangePolicy = Kokkos::RangePolicy<typename RowPtrType::execution_space>;

        const auto col_inds = IndicesType(
            Kokkos::view_alloc("col_inds", Kokkos::WithoutInitializing),
            static_cast<size_t>(num_non_zero)
        );

        const auto num_nodes = active_dofs.extent(0);
        const auto num_constraints = row_range.extent(0);

        Kokkos::parallel_for(
            "ComputeSystemColInds", RangePolicy(0, num_nodes),
            ComputeSystemColInds<RowPtrType, IndicesType>{
                num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
                node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
                target_node_freedom_table, row_range, row_ptrs, col_inds
            }
        );

        Kokkos::parallel_for(
            "ComputeConstraintsColInds", RangePolicy(0, num_constraints),
            ComputeConstraintsColInds<RowPtrType, IndicesType>{
                num_system_dofs, base_active_dofs, target_active_dofs, base_node_freedom_table,
                target_node_freedom_table, row_range, row_ptrs, col_inds
            }
        );

        return col_inds;
    }
};
}  // namespace kynema_fmb::solver
