#include <array>
#include <cstddef>
#include <ranges>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "solver/compute_row_ptrs.hpp"

namespace {

template <typename ValueType, typename DataType>
typename Kokkos::View<ValueType>::const_type CreateView(
    const std::string& name, const DataType& data
) {
    const auto view = Kokkos::View<ValueType>(Kokkos::view_alloc(name, Kokkos::WithoutInitializing));
    const auto host = typename Kokkos::View<ValueType, Kokkos::HostSpace>::const_type(data.data());
    const auto mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, view);
    Kokkos::deep_copy(mirror, host);
    Kokkos::deep_copy(view, mirror);
    return view;
}

}  // namespace

namespace kynema::solver::tests {

TEST(ComputeRowPtrs, OneElementOneNode) {
    constexpr auto num_system_dofs = 6U;
    constexpr auto num_dofs = num_system_dofs;

    const auto active_dofs = CreateView<size_t[1]>("active_dofs", std::array{6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[1]>("node_freedom_map_table", std::array{0UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nods_per_element", std::array{1UL});
    const auto node_state_indices = CreateView<size_t[1][1]>("node_state_indices", std::array{0UL});

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto row_ptrs = ComputeRowPtrs<Kokkos::View<size_t*>>::invoke(
        num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range
    );

    const auto row_ptrs_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_ptrs);

    for (auto row : std::views::iota(0U, 7U)) {
        EXPECT_EQ(row_ptrs_mirror(row), row * 6U);
    }
}

TEST(ComputeRowPtrs, OneElementTwoNodes) {
    constexpr auto num_system_dofs = 12U;
    constexpr auto num_dofs = num_system_dofs;

    const auto active_dofs = CreateView<size_t[2]>("active_dofs", std::array{6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[2]>("node_freedom_map_table", std::array{0UL, 6UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nods_per_element", std::array{2UL});
    const auto node_state_indices =
        CreateView<size_t[1][2]>("node_state_indices", std::array{0UL, 1UL});

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto row_ptrs = ComputeRowPtrs<Kokkos::View<size_t*>>::invoke(
        num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range
    );

    const auto row_ptrs_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_ptrs);

    for (auto row : std::views::iota(0U, 13U)) {
        EXPECT_EQ(row_ptrs_mirror(row), row * 12U);
    }
}

TEST(ComputeRowPtrs, TwoElementTwoNodesNoOverlap) {
    constexpr auto num_system_dofs = 24U;
    constexpr auto num_dofs = num_system_dofs;

    const auto active_dofs = CreateView<size_t[4]>("active_dofs", std::array{6UL, 6UL, 6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[4]>("node_freedom_map_table", std::array{0UL, 6UL, 12UL, 18UL});
    const auto num_nodes_per_element =
        CreateView<size_t[2]>("num_nods_per_element", std::array{2UL, 2UL});
    const auto node_state_indices =
        CreateView<size_t[2][2]>("node_state_indices", std::array{0UL, 1UL, 2UL, 3UL});

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto row_ptrs = ComputeRowPtrs<Kokkos::View<size_t*>>::invoke(
        num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range
    );

    const auto row_ptrs_mirror = Kokkos::create_mirror_view(row_ptrs);
    Kokkos::deep_copy(row_ptrs_mirror, row_ptrs);

    for (auto row : std::views::iota(0U, 25U)) {
        EXPECT_EQ(row_ptrs_mirror(row), row * 12U);
    }
}

TEST(ComputeRowPtrs, TwoElementTwoNodesOverlap) {
    constexpr auto num_system_dofs = 18U;
    constexpr auto num_dofs = num_system_dofs;

    const auto active_dofs = CreateView<size_t[3]>("active_dofs", std::array{6UL, 6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[3]>("node_freedom_map_table", std::array{0UL, 6UL, 12UL});
    const auto num_nodes_per_element =
        CreateView<size_t[2]>("num_nods_per_element", std::array{2UL, 2UL});
    const auto node_state_indices =
        CreateView<size_t[2][2]>("node_state_indices", std::array{0UL, 1UL, 1UL, 2UL});

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto row_ptrs = ComputeRowPtrs<Kokkos::View<size_t*>>::invoke(
        num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range
    );

    const auto row_ptrs_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_ptrs);

    for (auto row : std::views::iota(0U, 7U)) {
        EXPECT_EQ(row_ptrs_mirror(row), row * 12U);
    }

    for (auto row : std::views::iota(0U, 6U)) {
        EXPECT_EQ(row_ptrs_mirror(row + 7U), 72U + ((row + 1) * 18U));
    }

    for (auto row : std::views::iota(0U, 6U)) {
        EXPECT_EQ(row_ptrs_mirror(row + 13U), 180U + ((row + 1) * 12U));
    }
}

TEST(ComputeRowPtrs, OneElementOneNode_OneConstraint) {
    constexpr auto num_system_dofs = 12U;
    constexpr auto num_constraint_dofs = 6U;
    constexpr auto num_dofs = num_system_dofs + num_constraint_dofs;

    const auto active_dofs = CreateView<size_t[2]>("active_dofs", std::array{6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[2]>("node_freedom_map_table", std::array{0UL, 6UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nods_per_element", std::array{1UL});
    const auto node_state_indices = CreateView<size_t[1][1]>("node_state_indices", std::array{0UL});
    const auto base_active_dofs = CreateView<size_t[1]>("base_active_dofs", std::array{6UL});
    const auto target_active_dofs = CreateView<size_t[1]>("target_active_dofs", std::array{6UL});
    const auto base_node_freedom_table = CreateView<size_t[1][6]>(
        "base_node_freedom_table", std::array{0UL, 1UL, 2UL, 3UL, 4UL, 5UL}
    );
    const auto target_node_freedom_table = CreateView<size_t[1][6]>(
        "target_node_freedom_table", std::array{6UL, 7UL, 8UL, 9UL, 10UL, 11UL}
    );
    const auto row_range = CreateView<Kokkos::pair<size_t, size_t>[1]>(
        "row_range", std::array<Kokkos::pair<size_t, size_t>, 1>{Kokkos::make_pair(0UL, 6UL)}
    );

    const auto row_ptrs = ComputeRowPtrs<Kokkos::View<size_t*>>::invoke(
        num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range
    );

    const auto row_ptrs_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_ptrs);

    for (auto row : std::views::iota(0U, num_dofs + 1UL)) {
        EXPECT_EQ(row_ptrs_mirror(row), row * 12U);
    }
}

TEST(ComputeRowPtrs, OneElementOneNode_TwoConstraint) {
    constexpr auto num_system_dofs = 18U;
    constexpr auto num_constraint_dofs = 12U;
    constexpr auto num_dofs = num_system_dofs + num_constraint_dofs;

    const auto active_dofs = CreateView<size_t[3]>("active_dofs", std::array{6UL, 6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[3]>("node_freedom_map_table", std::array{0UL, 6UL, 12UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nods_per_element", std::array{1UL});
    const auto node_state_indices = CreateView<size_t[1][1]>("node_state_indices", std::array{0UL});
    const auto base_active_dofs = CreateView<size_t[2]>("base_active_dofs", std::array{6UL, 6UL});
    const auto target_active_dofs =
        CreateView<size_t[2]>("target_active_dofs", std::array{6UL, 6UL});
    const auto base_node_freedom_table = CreateView<size_t[2][6]>(
        "base_node_freedom_table",
        std::array{0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 0UL, 1UL, 2UL, 3UL, 4UL, 5UL}
    );
    const auto target_node_freedom_table = CreateView<size_t[2][6]>(
        "target_node_freedom_table",
        std::array{6UL, 7UL, 8UL, 9UL, 10UL, 11UL, 12UL, 13UL, 14UL, 15UL, 16UL, 17UL}
    );
    const auto row_range = CreateView<Kokkos::pair<size_t, size_t>[2]>(
        "row_range", std::array{Kokkos::make_pair(0UL, 6UL), Kokkos::make_pair(6UL, 12UL)}
    );

    const auto row_ptrs = ComputeRowPtrs<Kokkos::View<size_t*>>::invoke(
        num_system_dofs, num_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range
    );

    const auto row_ptrs_mirror = Kokkos::create_mirror_view(row_ptrs);
    Kokkos::deep_copy(row_ptrs_mirror, row_ptrs);

    for (auto row : std::views::iota(0U, 6U)) {
        EXPECT_EQ(row_ptrs_mirror(row), row * 18U);
    }

    for (auto row : std::views::iota(7U, num_dofs + 1UL)) {
        EXPECT_EQ(row_ptrs_mirror(row), (6U * 18U) + ((row - 6U) * 12U));
    }
}

}  // namespace kynema::solver::tests
