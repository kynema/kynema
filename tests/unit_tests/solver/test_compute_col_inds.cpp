#include <array>
#include <cstddef>
#include <ranges>
#include <string>

#include <KokkosSparse_SortCrs.hpp>
#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "solver/compute_col_inds.hpp"

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

TEST(ComputeColInds, OneElementOneNode) {
    constexpr auto num_non_zero = 36UL;
    constexpr auto num_system_dofs = 6UL;

    const auto active_dofs = CreateView<size_t[1]>("active_dofs", std::array{6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[1]>("node_freedom_map_table", std::array{0UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nodes_per_element", std::array{1UL});
    const auto node_state_indices = CreateView<size_t[1][1]>("node_state_indices", std::array{0UL});
    const auto row_ptrs = CreateView<int[7]>("row_ptrs", std::array{0, 6, 12, 18, 24, 30, 36});

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto col_inds = ComputeColInds<Kokkos::View<int[7]>, Kokkos::View<int*>>::invoke(
        num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range, row_ptrs
    );

    const auto col_inds_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), col_inds);

    for (auto row : std::views::iota(0U, 6U)) {
        for (auto col : std::views::iota(0U, 6U)) {
            EXPECT_EQ(col_inds_mirror((row * 6U) + col), col);
        }
    }
}

TEST(ComputeColInds, OneElementTwoNodes) {
    constexpr auto num_non_zero = 144UL;
    constexpr auto num_system_dofs = 12UL;

    const auto active_dofs = CreateView<size_t[2]>("active_dofs", std::array{6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[2]>("node_freedom_map_table", std::array{0UL, 6UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nodes_per_element", std::array{2UL});
    const auto node_state_indices =
        CreateView<size_t[1][2]>("node_state_indices", std::array{0UL, 1UL});
    const auto row_ptrs = CreateView<int[13]>(
        "row_ptrs", std::array{0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144}
    );

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto col_inds = ComputeColInds<Kokkos::View<int[13]>, Kokkos::View<int*>>::invoke(
        num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range, row_ptrs
    );

    KokkosSparse::sort_crs_graph(row_ptrs, col_inds);

    const auto col_inds_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), col_inds);

    for (auto row : std::views::iota(0U, 12U)) {
        for (auto col : std::views::iota(0U, 12U)) {
            EXPECT_EQ(col_inds_mirror((row * 12U) + col), col);
        }
    }
}

TEST(ComputeColInds, TwoElementsTwoNodesNoOverlap) {
    constexpr auto num_non_zero = 288UL;
    constexpr auto num_system_dofs = 24UL;

    const auto active_dofs = CreateView<size_t[4]>("active_dofs", std::array{6UL, 6UL, 6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[4]>("node_freedom_map_table", std::array{0UL, 6UL, 12UL, 18UL});
    const auto num_nodes_per_element =
        CreateView<size_t[2]>("num_nodes_per_element", std::array{2UL, 2UL});
    const auto node_state_indices =
        CreateView<size_t[2][2]>("node_state_indices", std::array{0UL, 1UL, 2UL, 3UL});
    const auto row_ptrs =
        CreateView<int[25]>("row_ptrs", std::array{0,   12,  24,  36,  48,  60,  72,  84,  96,
                                                   108, 120, 132, 144, 156, 168, 180, 192, 204,
                                                   216, 228, 240, 252, 264, 276, 288});

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto col_inds = ComputeColInds<Kokkos::View<int[25]>, Kokkos::View<int*>>::invoke(
        num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range, row_ptrs
    );

    KokkosSparse::sort_crs_graph(row_ptrs, col_inds);

    const auto col_inds_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), col_inds);

    for (auto row : std::views::iota(0U, 12U)) {
        for (auto col : std::views::iota(0U, 12U)) {
            EXPECT_EQ(col_inds_mirror((row * 12U) + col), col);
        }
    }

    for (auto row : std::views::iota(0U, 12U)) {
        for (auto col : std::views::iota(0U, 12U)) {
            EXPECT_EQ(col_inds_mirror(144U + (row * 12U) + col), 12U + col);
        }
    }
}

TEST(ComputeColInds, TwoElementsTwoNodesOverlap) {
    constexpr auto num_non_zero = 252UL;
    constexpr auto num_system_dofs = 18UL;

    const auto active_dofs = CreateView<size_t[3]>("active_dofs", std::array{6UL, 6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[3]>("node_freedom_map_table", std::array{0UL, 6UL, 12UL});
    const auto num_nodes_per_element =
        CreateView<size_t[2]>("num_nodes_per_element", std::array{2UL, 2UL});
    const auto node_state_indices =
        CreateView<size_t[2][2]>("node_state_indices", std::array{0UL, 1UL, 1UL, 2UL});
    const auto row_ptrs = CreateView<int[19]>(
        "row_ptrs",
        std::array{
            0, 12, 24, 36, 48, 60, 72, 90, 108, 126, 144, 162, 180, 192, 204, 216, 228, 240, 252
        }
    );

    const auto base_active_dofs = Kokkos::View<size_t*>("base_active_dofs", 0);
    const auto target_active_dofs = Kokkos::View<size_t*>("target_active_dofs", 0);
    const auto base_node_freedom_table = Kokkos::View<size_t* [6]>("base_node_freedom_table", 0);
    const auto target_node_freedom_table = Kokkos::View<size_t* [6]>("target_node_freedom_table", 0);
    const auto row_range = Kokkos::View<Kokkos::pair<size_t, size_t>*>("row_range", 0);

    const auto col_inds = ComputeColInds<Kokkos::View<int[19]>, Kokkos::View<int*>>::invoke(
        num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range, row_ptrs
    );

    KokkosSparse::sort_crs_graph(row_ptrs, col_inds);

    const auto col_inds_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), col_inds);

    for (auto row : std::views::iota(0U, 6U)) {
        for (auto col : std::views::iota(0U, 12U)) {
            EXPECT_EQ(col_inds_mirror((row * 12U) + col), col);
        }
    }

    for (auto row : std::views::iota(0U, 6U)) {
        for (auto col : std::views::iota(0U, 12U)) {
            EXPECT_EQ(col_inds_mirror(72U + (row * 18U) + col), col);
        }
    }

    for (auto row : std::views::iota(0U, 6U)) {
        for (auto col : std::views::iota(0U, 12U)) {
            EXPECT_EQ(col_inds_mirror(180U + (row * 12U) + col), 6U + col);
        }
    }
}

TEST(ComputeColInds, OneElementOneNode_OneConstraint) {
    constexpr auto num_non_zero = 216;
    constexpr auto num_system_dofs = 12U;

    const auto active_dofs = CreateView<size_t[2]>("active_dofs", std::array{6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[2]>("node_freedom_map_table", std::array{0UL, 6UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nodes_per_element", std::array{1UL});
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
    const auto row_ptrs = CreateView<int[19]>(
        "row_ptrs",
        std::array{
            0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216
        }
    );

    const auto col_inds = ComputeColInds<Kokkos::View<int[19]>, Kokkos::View<int*>>::invoke(
        num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range, row_ptrs
    );

    KokkosSparse::sort_crs_graph(row_ptrs, col_inds);

    const auto col_inds_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), col_inds);
    const auto row_ptrs_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_ptrs);

    for (auto row : std::views::iota(0, 6)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 0UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 1UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 2UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 3UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 4UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 5UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 12UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 13UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 14UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 15UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 16UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 17UL);
    }
    for (auto row : std::views::iota(6, 12)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 6UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 7UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 8UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 9UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 10UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 11UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 12UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 13UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 14UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 15UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 16UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 17UL);
    }
    for (auto row : std::views::iota(12, 18)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 0UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 1UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 2UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 3UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 4UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 5UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 6UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 7UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 8UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 9UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 10UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 11UL);
    }
}

TEST(ComputeColInds, OneElementOneNode_TwoConstraint) {
    constexpr auto num_non_zero = 396;
    constexpr auto num_system_dofs = 18U;

    const auto active_dofs = CreateView<size_t[3]>("active_dofs", std::array{6UL, 6UL, 6UL});
    const auto node_freedom_map_table =
        CreateView<size_t[3]>("node_freedom_map_table", std::array{0UL, 6UL, 12UL});
    const auto num_nodes_per_element =
        CreateView<size_t[1]>("num_nodes_per_element", std::array{1UL});
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
    const auto row_ptrs = CreateView<int[31]>(
        "row_ptrs",
        std::array{0,   18,  36,  54,  72,  90,  108, 120, 132, 144, 156, 168, 180, 192, 204, 216,
                   228, 240, 252, 264, 276, 288, 300, 312, 324, 336, 348, 360, 372, 384, 396}
    );

    const auto col_inds = ComputeColInds<Kokkos::View<int[31]>, Kokkos::View<int*>>::invoke(
        num_non_zero, num_system_dofs, active_dofs, node_freedom_map_table, num_nodes_per_element,
        node_state_indices, base_active_dofs, target_active_dofs, base_node_freedom_table,
        target_node_freedom_table, row_range, row_ptrs
    );

    KokkosSparse::sort_crs_graph(row_ptrs, col_inds);

    const auto col_inds_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), col_inds);
    const auto row_ptrs_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), row_ptrs);

    for (auto row : std::views::iota(0, 6)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 0UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 1UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 2UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 3UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 4UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 5UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 18UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 19UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 20UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 21UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 22UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 23UL);
        EXPECT_EQ(col_inds_mirror(entry + 12), 24UL);
        EXPECT_EQ(col_inds_mirror(entry + 13), 25UL);
        EXPECT_EQ(col_inds_mirror(entry + 14), 26UL);
        EXPECT_EQ(col_inds_mirror(entry + 15), 27UL);
        EXPECT_EQ(col_inds_mirror(entry + 16), 28UL);
        EXPECT_EQ(col_inds_mirror(entry + 17), 29UL);
    }
    for (auto row : std::views::iota(6, 12)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 6UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 7UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 8UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 9UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 10UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 11UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 18UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 19UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 20UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 21UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 22UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 23UL);
    }
    for (auto row : std::views::iota(12, 18)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 12UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 13UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 14UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 15UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 16UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 17UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 24UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 25UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 26UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 27UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 28UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 29UL);
    }
    for (auto row : std::views::iota(18, 24)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 0UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 1UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 2UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 3UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 4UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 5UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 6UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 7UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 8UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 9UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 10UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 11UL);
    }
    for (auto row : std::views::iota(24, 30)) {
        auto entry = row_ptrs_host(row);
        EXPECT_EQ(col_inds_mirror(entry), 0UL);
        EXPECT_EQ(col_inds_mirror(entry + 1), 1UL);
        EXPECT_EQ(col_inds_mirror(entry + 2), 2UL);
        EXPECT_EQ(col_inds_mirror(entry + 3), 3UL);
        EXPECT_EQ(col_inds_mirror(entry + 4), 4UL);
        EXPECT_EQ(col_inds_mirror(entry + 5), 5UL);
        EXPECT_EQ(col_inds_mirror(entry + 6), 12UL);
        EXPECT_EQ(col_inds_mirror(entry + 7), 13UL);
        EXPECT_EQ(col_inds_mirror(entry + 8), 14UL);
        EXPECT_EQ(col_inds_mirror(entry + 9), 15UL);
        EXPECT_EQ(col_inds_mirror(entry + 10), 16UL);
        EXPECT_EQ(col_inds_mirror(entry + 11), 17UL);
    }
}

}  // namespace kynema::solver::tests
