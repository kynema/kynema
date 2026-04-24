#pragma once

#include <KokkosSparse.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "compute_col_inds.hpp"
#include "compute_row_ptrs.hpp"

namespace kynema_fmb::solver {

/**
 * @brief The top level function object which creates the CRS matrix structure for the linear system
 * to be solved.
 *
 * @details This involves computing both he row pointers and the column indices for the entire
 * matrix. These indices will be unique and sorted
 */
template <typename CrsMatrixType>
struct CreateFullMatrix {
    using DeviceType = typename CrsMatrixType::device_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    [[nodiscard]] static CrsMatrixType invoke(
        size_t num_system_dofs, size_t num_dofs, const ConstView<size_t*>& base_active_dofs,
        const ConstView<size_t*>& target_active_dofs,
        const ConstView<size_t* [6]>& base_node_freedom_table,
        const ConstView<size_t* [6]>& target_node_freedom_table,
        const ConstView<Kokkos::pair<size_t, size_t>*>& row_range,
        const ConstView<size_t*>& active_dofs, const ConstView<size_t*>& node_freedom_map_table,
        const ConstView<size_t*>& num_nodes_per_element,
        const ConstView<size_t**>& node_state_indices
    ) {
        auto region = Kokkos::Profiling::ScopedRegion("Create Full Matrix");

        using ValuesType = typename CrsMatrixType::values_type::non_const_type;
        using RowPtrType = typename CrsMatrixType::staticcrsgraph_type::row_map_type::non_const_type;
        using IndicesType =
            typename CrsMatrixType::staticcrsgraph_type::entries_type::non_const_type;

        using RowPtrValueType = typename RowPtrType::value_type;
        using IndicesValueType = typename IndicesType::value_type;

        const auto row_ptrs = ComputeRowPtrs<RowPtrType>::invoke(
            num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
            node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
            target_node_freedom_table, row_range
        );

        auto num_non_zero = RowPtrValueType{};
        Kokkos::deep_copy(num_non_zero, Kokkos::subview(row_ptrs, num_dofs));

        const auto col_inds = ComputeColInds<RowPtrType, IndicesType>::invoke(
            num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table,
            num_nodes_per_element, node_state_indices, base_active_dofs, target_active_dofs,
            base_node_freedom_table, target_node_freedom_table, row_range, row_ptrs
        );

        KokkosSparse::sort_crs_graph(row_ptrs, col_inds);

        return CrsMatrixType(
            "full_matrix", static_cast<IndicesValueType>(num_dofs),
            static_cast<IndicesValueType>(num_dofs), num_non_zero,
            ValuesType(
                Kokkos::view_alloc("values", Kokkos::WithoutInitializing),
                static_cast<size_t>(num_non_zero)
            ),
            row_ptrs, col_inds
        );
    }
};
}  // namespace kynema_fmb::solver
