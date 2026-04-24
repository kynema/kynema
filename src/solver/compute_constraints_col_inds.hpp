#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel for computing the system constraints' contribution to the
 * column indicies for the CRS matrix to be solved at each nonlinear iteration
 */
template <typename RowPtrType, typename IndicesType>
struct ComputeConstraintsColInds {
    using RowPtrValueType = typename RowPtrType::value_type;
    using IndicesValueType = typename IndicesType::value_type;
    using DeviceType = typename RowPtrType::device_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t num_system_dofs{};
    ConstView<size_t*> base_active_dofs;
    ConstView<size_t*> target_active_dofs;
    ConstView<size_t* [6]> base_node_freedom_table;
    ConstView<size_t* [6]> target_node_freedom_table;
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    typename RowPtrType::const_type row_ptrs;
    IndicesType col_inds;

    KOKKOS_FUNCTION
    RowPtrValueType CalculateTargetInds(size_t constraint, RowPtrValueType dof_index) const {
        const auto target_cols = target_active_dofs(constraint);
        for (auto column = 0U; column < target_cols; ++column, ++dof_index) {
            const auto index = target_node_freedom_table(constraint, column);
            col_inds(dof_index) = static_cast<IndicesValueType>(index);
        }
        return dof_index;
    }

    KOKKOS_FUNCTION
    RowPtrValueType CalculateBaseInds(size_t constraint, RowPtrValueType dof_index) const {
        const auto base_cols = base_active_dofs(constraint);
        for (auto column = 0U; column < base_cols; ++column, ++dof_index) {
            const auto index = base_node_freedom_table(constraint, column);
            col_inds(dof_index) = static_cast<IndicesValueType>(index);
        }
        return dof_index;
    }

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        const auto first_row = row_range(constraint).first;
        const auto end_row = row_range(constraint).second;
        for (auto row = first_row; row < end_row; ++row) {
            const auto initial_index = row_ptrs(row + num_system_dofs);
            const auto second_index = CalculateTargetInds(constraint, initial_index);
            CalculateBaseInds(constraint, second_index);
        }
    }
};

}  // namespace kynema_fmb::solver
