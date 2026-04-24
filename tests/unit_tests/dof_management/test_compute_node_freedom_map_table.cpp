#include <array>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "dof_management/compute_node_freedom_map_table.hpp"
#include "dof_management/freedom_signature.hpp"
#include "state/state.hpp"

namespace kynema_fmb::dof::tests {

TEST(TestComputeNodeFreedomMapTable, OneNode) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(1U);
    Kokkos::deep_copy(state.node_freedom_allocation_table, FreedomSignature::AllComponents);

    compute_node_freedom_map_table(state);

    const auto host_node_freedom_map_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.node_freedom_map_table);

    EXPECT_EQ(host_node_freedom_map_table(0), 0);
}

TEST(TestComputeNodeFreedomMapTable, FourNodes) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(4U);
    constexpr auto host_node_freedom_allocation_table_data = std::array{
        FreedomSignature::AllComponents, FreedomSignature::JustPosition,
        FreedomSignature::JustRotation, FreedomSignature::NoComponents
    };
    const auto host_node_freedom_allocation_table =
        Kokkos::View<FreedomSignature[4], Kokkos::HostSpace>::const_type(
            host_node_freedom_allocation_table_data.data()
        );
    const auto mirror_node_freedom_allocation_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_allocation_table);
    Kokkos::deep_copy(mirror_node_freedom_allocation_table, host_node_freedom_allocation_table);
    Kokkos::deep_copy(state.node_freedom_allocation_table, mirror_node_freedom_allocation_table);

    compute_node_freedom_map_table(state);

    const auto host_node_freedom_map_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.node_freedom_map_table);

    EXPECT_EQ(host_node_freedom_map_table(0), 0);
    EXPECT_EQ(host_node_freedom_map_table(1), 6);
    EXPECT_EQ(host_node_freedom_map_table(2), 9);
    EXPECT_EQ(host_node_freedom_map_table(3), 12);
}

}  // namespace kynema_fmb::dof::tests
