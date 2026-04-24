#include <gtest/gtest.h>

#include "elements/elements.hpp"

namespace kynema_fmb::tests {

TEST(ElementsTest, DefaultConstructor) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    const auto elements = Elements<DeviceType>();
    EXPECT_EQ(elements.NumElementsInSystem(), 0);
    EXPECT_EQ(elements.beams.num_elems, 0);
    EXPECT_EQ(elements.masses.num_elems, 0);
    EXPECT_EQ(elements.springs.num_elems, 0);
}

TEST(ElementsTest, ConstructorWithBeams) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(1U, 2U, 2U);  // 1 beam element with 2 nodes, 2 qps
    auto masses = Masses<DeviceType>(0U);        // 0 mass elements
    auto springs = Springs<DeviceType>(0U);      // 0 spring elements
    const auto elements = Elements<DeviceType>(beams, masses, springs);
    EXPECT_EQ(elements.NumElementsInSystem(), 1);
    EXPECT_EQ(elements.beams.num_elems, 1);
    EXPECT_EQ(elements.masses.num_elems, 0);
    EXPECT_EQ(elements.springs.num_elems, 0);
}

TEST(ElementsTest, ConstructorWithMasses) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(0U, 0U, 0U);  // 0 beam elements
    auto masses = Masses<DeviceType>(1U);        // 1 mass element
    auto springs = Springs<DeviceType>(0U);      // 0 spring elements
    const auto elements = Elements<DeviceType>(beams, masses, springs);
    EXPECT_EQ(elements.NumElementsInSystem(), 1);
    EXPECT_EQ(elements.beams.num_elems, 0);
    EXPECT_EQ(elements.masses.num_elems, 1);
    EXPECT_EQ(elements.springs.num_elems, 0);
}

TEST(ElementsTest, ConstructorWithSprings) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(0U, 0U, 0U);  // 0 beam elements
    auto masses = Masses<DeviceType>(0U);        // 0 mass elements
    auto springs = Springs<DeviceType>(1U);      // 1 spring element
    const auto elements = Elements<DeviceType>(beams, masses, springs);
    EXPECT_EQ(elements.NumElementsInSystem(), 1);
    EXPECT_EQ(elements.beams.num_elems, 0);
    EXPECT_EQ(elements.masses.num_elems, 0);
    EXPECT_EQ(elements.springs.num_elems, 1);
}

TEST(ElementsTest, ConstructorWithBeamsMasses) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(1, 2, 2);  // 1 beam element with 2 nodes, 2 qps
    auto masses = Masses<DeviceType>(1);      // 1 mass element
    auto springs = Springs<DeviceType>(0U);   // 0 spring elements
    const auto elements = Elements<DeviceType>(beams, masses, springs);
    EXPECT_EQ(elements.NumElementsInSystem(), 2);
    EXPECT_EQ(elements.beams.num_elems, 1);
    EXPECT_EQ(elements.masses.num_elems, 1);
    EXPECT_EQ(elements.springs.num_elems, 0);
}

TEST(ElementsTest, ConstructorWithBeamsMassesSprings) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(1U, 2U, 2U);  // 1 beam element with 2 nodes, 2 qps
    auto masses = Masses<DeviceType>(1U);        // 1 mass element
    auto springs = Springs<DeviceType>(1U);      // 1 spring element
    const auto elements = Elements<DeviceType>(beams, masses, springs);
    EXPECT_EQ(elements.NumElementsInSystem(), 3U);
    EXPECT_EQ(elements.beams.num_elems, 1);
    EXPECT_EQ(elements.masses.num_elems, 1);
    EXPECT_EQ(elements.springs.num_elems, 1);
}

TEST(ElementsTest, NumberOfNodesPerElementBeams) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(2U, 3U, 2U);  // 2 beam elements with 3 nodes each
    Kokkos::deep_copy(beams.num_nodes_per_element, 3U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);
    const auto elements = Elements<DeviceType>{beams, masses, springs};

    EXPECT_EQ(elements.NumElementsInSystem(), 2);

    auto host_nodes_per_elem = Kokkos::create_mirror_view(elements.NumberOfNodesPerElement());
    Kokkos::deep_copy(host_nodes_per_elem, elements.NumberOfNodesPerElement());

    EXPECT_EQ(host_nodes_per_elem(0), 3);  // First beam element has 3 nodes
    EXPECT_EQ(host_nodes_per_elem(1), 3);  // Second beam element has 3 nodes
}

TEST(ElementsTest, NumberOfNodesPerElementMasses) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(3U);  // 3 mass elements
    auto springs = Springs<DeviceType>(0U);
    const auto elements = Elements<DeviceType>{beams, masses, springs};

    EXPECT_EQ(elements.NumElementsInSystem(), 3);

    auto host_nodes_per_elem = Kokkos::create_mirror_view(elements.NumberOfNodesPerElement());
    Kokkos::deep_copy(host_nodes_per_elem, elements.NumberOfNodesPerElement());

    // All mass elements should have 1 node
    EXPECT_EQ(host_nodes_per_elem(0), 1);
    EXPECT_EQ(host_nodes_per_elem(1), 1);
    EXPECT_EQ(host_nodes_per_elem(2), 1);
}

TEST(ElementsTest, NumberofNodesPerElementSprings) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(3U);  // 3 spring elements
    const auto elements = Elements<DeviceType>{beams, masses, springs};

    EXPECT_EQ(elements.NumElementsInSystem(), 3);

    auto host_nodes_per_elem = Kokkos::create_mirror_view(elements.NumberOfNodesPerElement());
    Kokkos::deep_copy(host_nodes_per_elem, elements.NumberOfNodesPerElement());

    // All spring elements should have 2 nodes
    EXPECT_EQ(host_nodes_per_elem(0), 2);
    EXPECT_EQ(host_nodes_per_elem(1), 2);
    EXPECT_EQ(host_nodes_per_elem(2), 2);
}

TEST(ElementsTest, NumberOfNodesPerElementBeamsMasses) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(2U, 4U, 2U);  // 2 beam elements with 4 nodes each
    Kokkos::deep_copy(beams.num_nodes_per_element, 4U);
    auto masses = Masses<DeviceType>(2U);    // 2 mass elements
    auto springs = Springs<DeviceType>(0U);  // 0 spring elements
    const auto elements = Elements<DeviceType>{beams, masses, springs};

    EXPECT_EQ(elements.NumElementsInSystem(), 4);

    auto host_nodes_per_elem = Kokkos::create_mirror_view(elements.NumberOfNodesPerElement());
    Kokkos::deep_copy(host_nodes_per_elem, elements.NumberOfNodesPerElement());

    // First two elements are beams
    EXPECT_EQ(host_nodes_per_elem(0), 4);
    EXPECT_EQ(host_nodes_per_elem(1), 4);

    // Last two elements are masses
    EXPECT_EQ(host_nodes_per_elem(2), 1);
    EXPECT_EQ(host_nodes_per_elem(3), 1);
}

TEST(ElementsTest, NodeStateIndicesBeams) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(2U, 3U, 2U);  // 2 beam elements with 3 nodes each
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(0U);

    // Set up state indices for beam nodes: element 1: [0,1,2], element 2: [2,3,4]
    auto host_beam_indices = Kokkos::create_mirror_view(beams.node_state_indices);
    host_beam_indices(0, 0) = 0;
    host_beam_indices(0, 1) = 1;
    host_beam_indices(0, 2) = 2;
    host_beam_indices(1, 0) = 2;
    host_beam_indices(1, 1) = 3;
    host_beam_indices(1, 2) = 4;
    Kokkos::deep_copy(beams.node_state_indices, host_beam_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 3U);

    const auto elements = Elements<DeviceType>{beams, masses, springs};
    auto indices = elements.NodeStateIndices();
    auto host_indices = Kokkos::create_mirror_view(indices);
    Kokkos::deep_copy(host_indices, indices);

    // Verify indices for both beam elements
    EXPECT_EQ(host_indices(0, 0), 0);
    EXPECT_EQ(host_indices(0, 1), 1);
    EXPECT_EQ(host_indices(0, 2), 2);
    EXPECT_EQ(host_indices(1, 0), 2);
    EXPECT_EQ(host_indices(1, 1), 3);
    EXPECT_EQ(host_indices(1, 2), 4);
}

TEST(ElementsTest, NodeStateIndicesMasses) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(3U);  // 3 mass elements
    auto springs = Springs<DeviceType>(0U);

    // Set up state indices for masses: [10, 20, 30]
    auto host_mass_indices = Kokkos::create_mirror_view(masses.state_indices);
    host_mass_indices(0) = 10;
    host_mass_indices(1) = 20;
    host_mass_indices(2) = 30;
    Kokkos::deep_copy(masses.state_indices, host_mass_indices);

    const auto elements = Elements<DeviceType>{beams, masses, springs};
    auto indices = elements.NodeStateIndices();
    auto host_indices = Kokkos::create_mirror_view(indices);
    Kokkos::deep_copy(host_indices, indices);

    // Verify indices for all mass elements (each mass has 1 node)
    EXPECT_EQ(host_indices(0, 0), 10);
    EXPECT_EQ(host_indices(1, 0), 20);
    EXPECT_EQ(host_indices(2, 0), 30);
}

TEST(ElementsTest, NodeStateIndicesSprings) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(0U, 0U, 0U);
    auto masses = Masses<DeviceType>(0U);
    auto springs = Springs<DeviceType>(1U);  // 1 spring element

    // Set up state indices for springs: [0,1]
    auto host_spring_indices = Kokkos::create_mirror_view(springs.node_state_indices);
    host_spring_indices(0, 0) = 0;
    host_spring_indices(0, 1) = 1;
    Kokkos::deep_copy(springs.node_state_indices, host_spring_indices);

    const auto elements = Elements<DeviceType>{beams, masses, springs};
    auto indices = elements.NodeStateIndices();
    auto host_indices = Kokkos::create_mirror_view(indices);
    Kokkos::deep_copy(host_indices, indices);
    EXPECT_EQ(host_indices(0, 0), 0);
    EXPECT_EQ(host_indices(0, 1), 1);
}

TEST(ElementsTest, NodeStateIndicesBeamsMassesSprings) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto beams = Beams<DeviceType>(1U, 2U, 2U);  // 1 beam element with 2 nodes
    auto masses = Masses<DeviceType>(2U);        // 2 mass elements
    auto springs = Springs<DeviceType>(1U);      // 1 spring element

    // Set up state indices for beam nodes: [0, 1]
    auto host_beam_indices = Kokkos::create_mirror_view(beams.node_state_indices);
    host_beam_indices(0, 0) = 0;
    host_beam_indices(0, 1) = 1;
    Kokkos::deep_copy(beams.node_state_indices, host_beam_indices);
    Kokkos::deep_copy(beams.num_nodes_per_element, 2U);

    // Set up state indices for masses: [10, 20]
    auto host_mass_indices = Kokkos::create_mirror_view(masses.state_indices);
    host_mass_indices(0) = 10;
    host_mass_indices(1) = 20;
    Kokkos::deep_copy(masses.state_indices, host_mass_indices);

    // Set up state indices for springs: [5, 15]
    auto host_spring_indices = Kokkos::create_mirror_view(springs.node_state_indices);
    host_spring_indices(0, 0) = 5;
    host_spring_indices(0, 1) = 15;
    Kokkos::deep_copy(springs.node_state_indices, host_spring_indices);

    const auto elements = Elements<DeviceType>{beams, masses, springs};
    auto indices = elements.NodeStateIndices();
    auto host_indices = Kokkos::create_mirror_view(indices);
    Kokkos::deep_copy(host_indices, indices);

    // Verify beam element indices - 1 element with 2 nodes
    EXPECT_EQ(host_indices(0, 0), 0);
    EXPECT_EQ(host_indices(0, 1), 1);

    // Verify mass element indices - 2 elements with 1 node each
    EXPECT_EQ(host_indices(1, 0), 10);
    EXPECT_EQ(host_indices(2, 0), 20);

    // Verify spring element indices - 1 element with 2 nodes
    EXPECT_EQ(host_indices(3, 0), 5);
    EXPECT_EQ(host_indices(3, 1), 15);
}

}  // namespace kynema_fmb::tests
