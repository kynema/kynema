#include <gtest/gtest.h>

#include "model/model.hpp"
#include "test_utilities.hpp"

namespace {
inline auto SetUpSprings() {
    auto model = kynema::Model();

    // Add two nodes for the spring element
    auto node1_id =
        model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();  // First node at origin
    auto node2_id =
        model.AddNode().SetPosition(1., 0., 0., 1., 0., 0., 0.).Build();  // Second node at (1,0,0)

    const auto k = 10.;  // Spring stiffness coeff.
    const auto l0 = 1.;  // Undeformed length
    model.AddSpringElement(node1_id, node2_id, k, l0);

    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    auto springs = model.CreateSprings<DeviceType>();
    auto state = model.CreateState<DeviceType>();

    return std::make_tuple(springs, state);
}

}  // namespace

namespace kynema::tests {

TEST(SpringsTest, NodeStateIndices) {
    auto [springs, _] = SetUpSprings();
    auto node_state_indices_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), springs.node_state_indices);
    EXPECT_EQ(node_state_indices_host(0, 0), 0);
    EXPECT_EQ(node_state_indices_host(0, 1), 1);
}

TEST(SpringsTest, InitialPositionVector) {
    auto [springs, _] = SetUpSprings();
    expect_kokkos_view_2D_equal(
        springs.x0,
        {
            {1., 0., 0.},  // Vector from node 0 to node 1
        }
    );
}

TEST(SpringsTest, ReferenceLength) {
    auto [springs, _] = SetUpSprings();
    expect_kokkos_view_1D_equal(springs.l_ref, {1.});  // Undeformed length
}

TEST(SpringsTest, SpringStiffness) {
    auto [springs, _] = SetUpSprings();
    expect_kokkos_view_1D_equal(springs.k, {10.});  // Spring stiffness
}

}  // namespace kynema::tests
