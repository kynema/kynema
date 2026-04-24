#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief Kernel to compute the constraints' contribution to the row pointers of the CRS matrix
 */
template <typename RowPtrType>
struct ComputeConstraintsRowEntries {
    using ValueType = typename RowPtrType::value_type;
    using DeviceType = typename RowPtrType::device_type;
    template <typename value_type>
    using View = Kokkos::View<value_type, DeviceType>;
    template <typename value_type>
    using ConstView = typename View<value_type>::const_type;

    size_t num_system_dofs{};
    ConstView<size_t*> base_active_dofs;
    ConstView<size_t*> target_active_dofs;
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    RowPtrType row_entries;

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        const auto first_row = row_range(constraint).first;
        const auto end_row = row_range(constraint).second;
        const auto num_base = static_cast<ValueType>(base_active_dofs(constraint));
        const auto num_target = static_cast<ValueType>(target_active_dofs(constraint));
        for (auto row = first_row; row < end_row; ++row) {
            row_entries(row + num_system_dofs) = num_base + num_target;
        }
    }
};

}  // namespace kynema_fmb::solver
