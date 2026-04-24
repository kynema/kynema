#include <array>
#include <cstddef>
#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "dof_management/create_element_freedom_table.hpp"
#include "elements/beams/beams.hpp"
#include "elements/elements.hpp"
#include "elements/masses/masses.hpp"
#include "elements/springs/springs.hpp"
#include "state/state.hpp"

namespace kynema_fmb::dof::tests {

TEST(TestCreateElementFreedomTable, OneBeamElementWithOneNode_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(1U);
    Kokkos::deep_copy(state.node_freedom_map_table, 0U);

    auto beams = Beams<DeviceType>(1U, 1U, 1U);
    Kokkos::deep_copy(beams.node_state_indices, 0U);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 1 beam element + 0 mass elements + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 1);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);

    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 0, k), k
        );  // Beam Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
}

TEST(TestCreateElementFreedomTable, OneMassElementWithOneNode_NoBeamNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(1U);
    Kokkos::deep_copy(state.node_freedom_map_table, 0U);

    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(1U);
    Kokkos::deep_copy(masses.state_indices, 0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 0 beam elements + 1 mass element + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 1);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), masses.element_freedom_table);

    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_element_freedom_table(0, k), k
        );  // Mass Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
}

TEST(TestCreateElementFreedomTable, OneSpringElementWithTwoNodes_NoBeamNoMass) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 3UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(1U);
    constexpr auto host_node_state_indices_data = std::array{0UL, 1UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[1][2], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, springs.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(springs.node_state_indices, mirror_node_state_indices);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 0 beam elements + 0 mass elements + 1 spring element

    EXPECT_EQ(elements.NumElementsInSystem(), 1);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), springs.element_freedom_table);

    for (auto k : std::views::iota(0, 3)) {
        EXPECT_EQ(host_element_freedom_table(0, 0, k), k);  // Spring Element 1 Nodes 1 DOFs: 0, 1, 2
    }
    for (auto k : std::views::iota(0U, 3U)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 1, k), k + 3U
        );  // Spring Element 1 Nodes 2 DOFs: 3, 4, 5
    }
}

TEST(TestCreateElementFreedomTable, OneBeamElementWithTwoNodes_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(1U, 2U, 1U);
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
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 1 beam element + 0 mass elements + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 1);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);

    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 0, k), k
        );  // Beam Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 1, k), k + 6U
        );  // Beam Element 1 Node 2 DOFs: 6, 7, 8, 9, 10, 11
    }
}

TEST(TestCreateElementFreedomTable, OneBeamElementWithOneNode_OneMassElementWithOneNode_NoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(1U, 1U, 1U);
    Kokkos::deep_copy(beams.node_state_indices, 0U);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);
    auto masses = Masses<DeviceType>(1U);
    Kokkos::deep_copy(masses.state_indices, 1U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 1 beam element + 1 mass element + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 2);

    create_element_freedom_table(elements, state);

    const auto host_beams_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_beams_element_freedom_table(0, 0, k), k
        );  // Beam Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }

    const auto host_masses_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), masses.element_freedom_table);
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_masses_element_freedom_table(0, k), k + 6U
        );  // Mass Element 1 Node 1 DOFs: 6, 7, 8, 9, 10, 11
    }
}

TEST(TestCreateElementFreedomTable, TwoBeamElementsWithOneNode_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(2U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(2U, 1U, 1U);
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
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 2 beam elements + 0 mass elements + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 2);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);

    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 0, k), k
        );  // Beam Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(1, 0, k), k + 6U
        );  // Beam Element 2 Node 1 DOFs: 6, 7, 8, 9, 10, 11
    }
}

TEST(TestCreateElementFreedomTable, TwoBeamElementsWithTwoNodesShared_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(3U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL, 12UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[3], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(2U, 2U, 1U);
    constexpr auto host_node_state_indices_data = std::array{0UL, 1UL, 1UL, 2UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[2][2], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, beams.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(beams.node_state_indices, mirror_node_state_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 2U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 2 beam elements + 0 mass elements + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 2);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);

    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 0, k), k
        );  // Beam Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 1, k), k + 6U
        );  // Beam Element 1 Node 2 DOFs: 6, 7, 8, 9, 10, 11
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(1, 0, k), k + 6U
        );  // Beam Element 2 Node 1 DOFs: 6, 7, 8, 9, 10, 11
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(1, 1, k), k + 12U
        );  // Beam Element 2 Node 2 DOFs: 12, 13, 14, 15, 16, 17
    }
}

TEST(TestCreateElementFreedomTable, TwoBeamElementsWithTwoNodesShared_Flipped_NoMassNoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(3U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL, 12UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[3], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(2U, 2U, 1U);
    constexpr auto host_node_state_indices_data = std::array{1UL, 2UL, 0UL, 1UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[2][2], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, beams.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(beams.node_state_indices, mirror_node_state_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 2U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 2 beam elements + 0 mass elements + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 2);

    create_element_freedom_table(elements, state);

    const auto host_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);

    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_element_freedom_table(1, 0, k), k
        );  // Beam Element 2 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(1, 1, k), k + 6U
        );  // Beam Element 2 Node 2 DOFs: 6, 7, 8, 9, 10, 11
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 0, k), k + 6U
        );  // Beam Element 1 Node 1 DOFs: 6, 7, 8, 9, 10, 11
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_element_freedom_table(0, 1, k), k + 12U
        );  // Beam Element 1 Node 2 DOFs: 12, 13, 14, 15, 16, 17
    }
}

TEST(TestCreateElementFreedomTable, TwoBeamElementsWithOneNode_OneMassElementWithOneNode_NoSpring) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(3U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL, 12UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[3], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(2U, 1U, 1U);
    constexpr auto host_node_state_indices_data = std::array{0UL, 1UL};
    const auto host_node_state_indices =
        Kokkos::View<size_t[2][1], Kokkos::HostSpace>::const_type(host_node_state_indices_data.data()
        );
    const auto mirror_node_state_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, beams.node_state_indices);
    Kokkos::deep_copy(mirror_node_state_indices, host_node_state_indices);
    Kokkos::deep_copy(beams.node_state_indices, mirror_node_state_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);
    auto masses = Masses<DeviceType>(1U);
    Kokkos::deep_copy(masses.state_indices, 2U);
    auto springs = Springs<DeviceType>(0U);
    auto elements = Elements<DeviceType>{
        beams, masses, springs
    };  // 2 beam elements + 1 mass element + 0 spring elements

    EXPECT_EQ(elements.NumElementsInSystem(), 3);

    create_element_freedom_table(elements, state);

    const auto host_beams_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_beams_element_freedom_table(0, 0, k), k
        );  // Beam Element 1 Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_beams_element_freedom_table(1, 0, k), k + 6U
        );  // Beam Element 2 Node 1 DOFs: 6, 7, 8, 9, 10, 11
    }

    const auto host_masses_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), masses.element_freedom_table);
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_masses_element_freedom_table(0, k), k + 12U
        );  // Mass Element 1 Node 1 DOFs: 12, 13, 14, 15, 16, 17
    }
}

TEST(
    TestCreateElementFreedomTable,
    OneBeamElementWithOneNode_OneMassElementWithOneNode_OneSpringElementWithTwoNodes
) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state = State<DeviceType>(4U);
    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL, 12UL, 15UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[4], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    auto beams = Beams<DeviceType>(1U, 1U, 1U);
    Kokkos::deep_copy(beams.node_state_indices, 0U);
    Kokkos::deep_copy(beams.num_nodes_per_element, 1U);

    auto masses = Masses<DeviceType>(1U);
    Kokkos::deep_copy(masses.state_indices, 1U);

    auto springs = Springs<DeviceType>(1U);
    constexpr auto host_spring_node_indices_data = std::array{2UL, 3UL};
    const auto host_spring_node_indices = Kokkos::View<size_t[1][2], Kokkos::HostSpace>::const_type(
        host_spring_node_indices_data.data()
    );
    const auto mirror_spring_node_indices =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, springs.node_state_indices);
    Kokkos::deep_copy(mirror_spring_node_indices, host_spring_node_indices);
    Kokkos::deep_copy(springs.node_state_indices, mirror_spring_node_indices);

    auto elements =
        Elements<DeviceType>{beams, masses, springs};  // 1 beam + 1 mass + 1 spring element

    EXPECT_EQ(elements.NumElementsInSystem(), 3);

    create_element_freedom_table(elements, state);

    const auto host_beams_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), beams.element_freedom_table);
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_beams_element_freedom_table(0, 0, k), k
        );  // Beam Element Node 1 DOFs: 0, 1, 2, 3, 4, 5
    }

    const auto host_masses_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), masses.element_freedom_table);
    for (auto k : std::views::iota(0U, 6U)) {
        EXPECT_EQ(
            host_masses_element_freedom_table(0, k), k + 6U
        );  // Mass Element Node 1 DOFs: 6, 7, 8, 9, 10, 11
    }

    const auto host_springs_element_freedom_table =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), springs.element_freedom_table);
    for (auto k : std::views::iota(0U, 3U)) {
        EXPECT_EQ(
            host_springs_element_freedom_table(0, 0, k), k + 12U
        );  // Spring Element Node 1 DOFs: 12, 13, 14
    }
    for (auto k : std::views::iota(0U, 3U)) {
        EXPECT_EQ(
            host_springs_element_freedom_table(0, 1, k), k + 15U
        );  // Spring Element Node 2 DOFs: 15, 16, 17
    }
}

}  // namespace kynema_fmb::dof::tests
