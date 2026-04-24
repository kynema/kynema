#include <gtest/gtest.h>

#include "dof_management/assemble_node_freedom_allocation_table.hpp"

namespace kynema_fmb::dof::tests {

TEST(TestAssembleNodeFreedomAllocationTable, OneBeamElementWithOneNode_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(1U);  // 1 node in the system

    auto beams = Beams<DeviceType>(1U, 1U, 1U);  // 1 beam element with 1 node per element
    Kokkos::deep_copy(beams.node_state_indices, 0U);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::AllComponents);
}

TEST(TestAssembleNodeFreedomAllocationTable, OneMassElementWithOneNode_NoBeamNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(1U);  // 1 node in the system

    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(1U);  // 1 mass element with 1 node
    Kokkos::deep_copy(masses.state_indices, 0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::AllComponents);
}

TEST(TestAssembleNodeFreedomAllocationTable, OneSpringElementWithTwoNodes_NoBeamNoMass) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);  // 2 nodes in the system

    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(1U);  // 1 spring element with 2 nodes
    constexpr auto host_node_state_indices_data = std::array{0UL, 1UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[1][2], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, springs.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(springs.node_state_indices, mirror_node_state_indices);
    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::JustPosition);
    EXPECT_EQ(host_node_freedom_allocation_table(1), FreedomSignature::JustPosition);
}

TEST(
    TestAssembleNodeFreedomAllocationTable,
    OneBeamElementWithOneNode_OneMassElementWithOneNode_NoSpring
) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);  // 2 nodes in the system

    auto beams = Beams<DeviceType>(1U, 1U, 1U);  // 1 beam element with 1 node per element
    Kokkos::deep_copy(beams.node_state_indices, 0U);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);
    auto masses = Masses<DeviceType>(1U);  // 1 mass element with 1 node
    Kokkos::deep_copy(masses.state_indices, 1U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::AllComponents);
    EXPECT_EQ(host_node_freedom_allocation_table(1), FreedomSignature::AllComponents);
}

TEST(TestAssembleNodeFreedomAllocationTable, OneBeamElementWithTwoNodes_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);  // 2 nodes in the system

    auto beams = Beams<DeviceType>(1U, 2U, 1U);  // 1 beam element with 2 nodes per element
    constexpr auto host_node_state_indices_data = std::array{0UL, 1UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[1][2], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, beams.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(beams.node_state_indices, mirror_node_state_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 2U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::AllComponents);
    EXPECT_EQ(host_node_freedom_allocation_table(1), FreedomSignature::AllComponents);
}

TEST(TestAssembleNodeFreedomAllocationTable, TwoBeamElementsWithOneNode_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);  // 2 nodes in the system

    auto beams = Beams<DeviceType>(2U, 1U, 1U);  // 2 beam elements with 1 node per element
    constexpr auto host_node_state_indices_data = std::array{0UL, 1UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[2][1], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, beams.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(beams.node_state_indices, mirror_node_state_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::AllComponents);
    EXPECT_EQ(host_node_freedom_allocation_table(1), FreedomSignature::AllComponents);
}

TEST(
    TestAssembleNodeFreedomAllocationTable,
    OneBeamElementWithTwoNodes_OneMassElementWithOneNode_OneSpringElementWithTwoNodes
) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(5U);  // 5 nodes in the system

    auto beams = Beams<DeviceType>(1U, 2U, 1U);  // 1 beam element with 2 nodes per element
    constexpr auto host_node_state_indices_data_beams = std::array{0UL, 1UL};
    const auto host_node_state_indices_beams =
        Kokkos::View<size_t[1][2], Kokkos::HostSpace>::const_type(
            host_node_state_indices_data_beams.data()
        );
    const auto mirror_node_state_indices_beams =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, beams.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices_beams, host_node_state_indices_beams);
    Kokkos::deep_copy(beams.node_state_indices, mirror_node_state_indices_beams);
    Kokkos::deep_copy(beams.num_nodes_per_element, 2U);

    auto masses = Masses<DeviceType>(1U);  // 1 mass element with 1 node
    Kokkos::deep_copy(masses.state_indices, 2U);

    auto springs = Springs<DeviceType>(1U);  // 1 spring element with 2 nodes
    constexpr auto host_node_state_indices_data_springs = std::array{3UL, 4UL};
    const auto host_node_state_indices_springs =
        Kokkos::View<size_t[1][2], Kokkos::HostSpace>::const_type(
            host_node_state_indices_data_springs.data()
        );
    const auto mirror_node_state_indices_springs =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, springs.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices_springs, host_node_state_indices_springs);
    Kokkos::deep_copy(springs.node_state_indices, mirror_node_state_indices_springs);

    auto elements = Elements<DeviceType>{beams, masses, springs};

    auto constraints =
        Constraints<DeviceType>(std::vector<constraints::Constraint>{}, std::vector<Node>{});

    assemble_node_freedom_allocation_table(state, elements, constraints);

    const auto host_node_freedom_allocation_table = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), state.node_freedom_allocation_table
    );

    EXPECT_EQ(host_node_freedom_allocation_table(0), FreedomSignature::AllComponents);
    EXPECT_EQ(host_node_freedom_allocation_table(1), FreedomSignature::AllComponents);
    EXPECT_EQ(host_node_freedom_allocation_table(2), FreedomSignature::AllComponents);
    EXPECT_EQ(host_node_freedom_allocation_table(3), FreedomSignature::JustPosition);
    EXPECT_EQ(host_node_freedom_allocation_table(4), FreedomSignature::JustPosition);
}

}  // namespace kynema_fmb::dof::tests
