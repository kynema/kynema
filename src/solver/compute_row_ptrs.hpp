#pragma once

#include <Kokkos_Core.hpp>

#include "compute_constraints_row_entries.hpp"
#include "compute_system_row_entries.hpp"
#include "scan_row_entries.hpp"

namespace kynema_fmb::solver {

/**
 * @brief Top level function object for calculating the row pointers of the CRS matrix to be solved
 * during each nonlinear iteration.
 */
template <typename RowPtrType>
struct ComputeRowPtrs {
    using DeviceType = typename RowPtrType::device_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    [[nodiscard]] static RowPtrType invoke(
        size_t num_system_dofs, size_t num_dofs, const ConstView<size_t*>& active_dofs,
        const ConstView<size_t*>& node_freedom_map_table,
        const ConstView<size_t*>& num_nodes_per_element,
        const ConstView<size_t**>& node_state_indices, const ConstView<size_t*>& base_active_dofs,
        const ConstView<size_t*>& target_active_dofs,
        const ConstView<size_t* [6]>& base_node_freedom_table,
        const ConstView<size_t* [6]>& target_node_freedom_table,
        const ConstView<Kokkos::pair<size_t, size_t>*>& row_range
    ) {
        using Kokkos::deep_copy;
        using Kokkos::parallel_for;
        using Kokkos::subview;
        using Kokkos::view_alloc;
        using Kokkos::WithoutInitializing;
        using RangePolicy = Kokkos::RangePolicy<typename RowPtrType::execution_space>;

        const auto row_entries =
            RowPtrType(view_alloc("row_entries", WithoutInitializing), num_dofs);
        const auto row_ptrs = RowPtrType(view_alloc("row_ptrs", WithoutInitializing), num_dofs + 1);
        deep_copy(subview(row_ptrs, 0), 0UL);

        const auto num_nodes = active_dofs.extent(0);
        const auto num_constraints = row_range.extent(0);

        parallel_for(
            "ComputeSystemRowEntries", RangePolicy(0, num_nodes),
            ComputeSystemRowEntries<RowPtrType>{
                active_dofs, node_freedom_map_table, num_nodes_per_element, node_state_indices,
                base_active_dofs, target_active_dofs, base_node_freedom_table,
                target_node_freedom_table, row_range, row_entries
            }
        );

        parallel_for(
            "ComputeConstraintsRowEntries", RangePolicy(0, num_constraints),
            ComputeConstraintsRowEntries<RowPtrType>{
                num_system_dofs, base_active_dofs, target_active_dofs, row_range, row_entries
            }
        );

        parallel_scan(
            "ScanRowEntries", RangePolicy(0, num_dofs),
            ScanRowEntries<RowPtrType>{row_entries, row_ptrs}
        );

        return row_ptrs;
    }
};
}  // namespace kynema_fmb::solver
