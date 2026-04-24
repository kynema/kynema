#pragma once

#include <KokkosSparse.hpp>
#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which copies the transpose gradient matrix contributions for the
 * base and target node of a constraint to the correct location in the global CRS matrix.
 */
template <typename CrsMatrixType>
struct CopyConstraintsTransposeToSparseMatrix {
    using DeviceType = typename CrsMatrixType::device_type;
    using RowDataType = typename CrsMatrixType::values_type::non_const_type;
    using ColIdxType = typename CrsMatrixType::staticcrsgraph_type::entries_type::non_const_type;
    using TeamPolicy = typename Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t num_system_cols{};
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
        constexpr bool is_sorted = true;
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        const auto num_cols =
            static_cast<int>(row_range(constraint).second - row_range(constraint).first);
        const auto first_col = static_cast<typename CrsMatrixType::ordinal_type>(
            row_range(constraint).first + num_system_cols
        );
        const auto num_base_dofs =
            static_cast<int>(count_active_dofs(base_node_freedom_signature(constraint)));
        const auto base_start_row = static_cast<int>(base_node_freedom_table(constraint, 0));
        const auto base_end_row = base_start_row + num_base_dofs;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, base_start_row, base_end_row),
            [&](int i) {
                const auto row_number = i - base_start_row;
                const auto hint = static_cast<typename CrsMatrixType::ordinal_type>(0);
                auto row = sparse.row(i);

                auto offset = KokkosSparse::findRelOffset(
                    &(row.colidx(0)), row.length, first_col, hint, is_sorted
                );
                for (auto entry = 0; entry < num_cols; ++entry, ++offset) {
                    if constexpr (force_atomic) {
                        Kokkos::atomic_add(
                            &row.value(offset), base_dense(constraint, row_number, entry)
                        );
                    } else {
                        row.value(offset) = base_dense(constraint, row_number, entry);
                    }
                }
            }
        );

        const auto num_target_dofs =
            static_cast<int>(count_active_dofs(target_node_freedom_signature(constraint)));
        const auto target_start_row = static_cast<int>(target_node_freedom_table(constraint, 0));
        const auto target_end_row = target_start_row + num_target_dofs;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, target_start_row, target_end_row),
            [&](int i) {
                const auto row_number = i - target_start_row;
                const auto hint = static_cast<typename CrsMatrixType::ordinal_type>(0);
                auto row = sparse.row(i);

                auto offset = KokkosSparse::findRelOffset(
                    &(row.colidx(0)), row.length, first_col, hint, is_sorted
                );
                for (auto entry = 0; entry < num_cols; ++entry, ++offset) {
                    if constexpr (force_atomic) {
                        Kokkos::atomic_add(
                            &row.value(offset), target_dense(constraint, row_number, entry)
                        );
                    } else {
                        row.value(offset) = target_dense(constraint, row_number, entry);
                    }
                }
            }
        );
    }
};
}  // namespace kynema_fmb::solver
