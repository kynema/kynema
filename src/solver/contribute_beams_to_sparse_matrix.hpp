#pragma once

#include <KokkosSparse.hpp>
#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which sums the system matrix contributions computed at each node in a beam
 * into the correct location of the global CRS matrix.
 */
template <typename CrsMatrixType>
struct ContributeBeamsToSparseMatrix {
    using DeviceType = typename CrsMatrixType::device_type;
    using RowDataType = typename CrsMatrixType::values_type::non_const_type;
    using ColIdxType = typename CrsMatrixType::staticcrsgraph_type::entries_type::non_const_type;
    using TeamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    double conditioner{};
    ConstView<size_t*> num_nodes_per_element;
    ConstView<dof::FreedomSignature**> element_freedom_signature;
    ConstView<size_t** [6]> element_freedom_table;
    ConstView<double*** [6][6]> dense;
    CrsMatrixType sparse;

    KOKKOS_FUNCTION
    void operator()(member_type member) const {
        const auto element = member.league_rank();
        const auto num_nodes = static_cast<int>(num_nodes_per_element(element));
        constexpr auto is_sorted = true;
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, num_nodes * num_nodes), [&](int node_12) {
                const auto node_1 = node_12 % num_nodes;
                const auto node_2 = node_12 / num_nodes;
                constexpr auto num_dofs = 6;
                const auto hint =
                    static_cast<typename CrsMatrixType::ordinal_type>(node_2 * num_dofs);

                const auto first_column = static_cast<typename CrsMatrixType::ordinal_type>(
                    element_freedom_table(element, node_2, 0)
                );
                const auto local_dense =
                    Kokkos::subview(dense, element, node_1, node_2, Kokkos::ALL, Kokkos::ALL);
                for (auto component_1 = 0; component_1 < num_dofs; ++component_1) {
                    const auto row_num =
                        static_cast<int>(element_freedom_table(element, node_1, component_1));
                    auto row = sparse.row(row_num);
                    auto offset = KokkosSparse::findRelOffset(
                        &(row.colidx(0)), row.length, first_column, hint, is_sorted
                    );
                    for (auto component_2 = 0; component_2 < num_dofs; ++component_2, ++offset) {
                        const auto contribution =
                            local_dense(component_1, component_2) * conditioner;
                        if constexpr (force_atomic) {
                            Kokkos::atomic_add(&(row.value(offset)), contribution);
                        } else {
                            row.value(offset) += contribution;
                        }
                    }
                }
            }
        );
    }
};

}  // namespace kynema_fmb::solver
