#pragma once

#include <KokkosSparse.hpp>
#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which copies the gradient matrix contributions for the base and target node of
 * a constraint to the correct locations in the global CRS matrix.
 */
template <typename CrsMatrixType>
struct CopyConstraintsToSparseMatrix {
    using DeviceType = typename CrsMatrixType::device_type;
    using RowDataType = typename CrsMatrixType::values_type::non_const_type;
    using ColIdxType = typename CrsMatrixType::staticcrsgraph_type::entries_type::non_const_type;
    using TeamPolicy = typename Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t num_system_rows{};
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    ConstView<dof::FreedomSignature*> base_node_freedom_signature;
    ConstView<dof::FreedomSignature*> target_node_freedom_signature;
    ConstView<size_t* [6]> base_node_freedom_table;
    ConstView<size_t* [6]> target_node_freedom_table;
    ConstView<double* [6][6]> base_dense;
    ConstView<double* [6][6]> target_dense;
    CrsMatrixType sparse;

    KOKKOS_FUNCTION
    void operator()(member_type member) const {
        const auto constraint = member.league_rank();
        constexpr auto is_sorted = true;
        const auto start_row = row_range(constraint).first;
        const auto end_row = row_range(constraint).second;
        const auto num_base_dofs = count_active_dofs(base_node_freedom_signature(constraint));
        const auto num_target_dofs = count_active_dofs(target_node_freedom_signature(constraint));
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;

        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, start_row, end_row), [&](size_t i) {
            const auto row_number = i - start_row;
            const auto real_row = num_system_rows + i;
            constexpr auto hint = static_cast<typename CrsMatrixType::ordinal_type>(0);
            auto row = sparse.row(static_cast<int>(real_row));
            auto first_col = static_cast<typename CrsMatrixType::ordinal_type>(
                base_node_freedom_table(constraint, 0)
            );
            auto offset = KokkosSparse::findRelOffset(
                &(row.colidx(0)), row.length, first_col, hint, is_sorted
            );
            for (auto entry = 0U; entry < num_base_dofs; ++entry, ++offset) {
                if constexpr (force_atomic) {
                    Kokkos::atomic_add(
                        &row.value(offset), base_dense(constraint, row_number, entry)
                    );
                } else {
                    row.value(offset) = base_dense(constraint, row_number, entry);
                }
            }

            first_col = static_cast<typename CrsMatrixType::ordinal_type>(
                target_node_freedom_table(constraint, 0)
            );
            offset = KokkosSparse::findRelOffset(
                &(row.colidx(0)), row.length, first_col, hint, is_sorted
            );
            for (auto entry = 0U; entry < num_target_dofs; ++entry, ++offset) {
                if constexpr (force_atomic) {
                    Kokkos::atomic_add(
                        &row.value(offset), target_dense(constraint, row_number, entry)
                    );
                } else {
                    row.value(offset) = target_dense(constraint, row_number, entry);
                }
            }
        });
    }
};
}  // namespace kynema_fmb::solver
