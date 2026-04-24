#include <gtest/gtest.h>

#include "dof_management/assemble_node_freedom_allocation_table.hpp"
#include "dof_management/create_constraint_freedom_table.hpp"

namespace kynema_fmb::dof::tests {

TEST(TestCreateConstraintFreedomTable, SingleNodeConstraint_FixedBC) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto invalid_node = Node(0U, {0., 0., 0., 1., 0., 0., 0.});  // base node - index is 0
    auto node_1 = Node(1U, {1., 0., 0., 1., 0., 0., 0.});        // target node - index is 1
    auto fixed_bc = constraints::Constraint(0, constraints::ConstraintType::FixedBC, {0, 1});
    auto constraints =
        Constraints<DeviceType>(std::array{fixed_bc}, std::array{invalid_node, node_1});

    auto elements = Elements<DeviceType>();
    auto state = State<DeviceType>(2U);  // 2 nodes in the system
    assemble_node_freedom_allocation_table(state, elements, constraints);

    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 17UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    create_constraint_freedom_table(constraints, state);

    const auto host_target_node_index =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_index);
    const auto host_target_node_freedom_signature =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_freedom_signature);
    const auto host_target_node_freedom_table =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_freedom_table);

    EXPECT_EQ(host_target_node_index(0), 1);  // target node index is 1
    EXPECT_EQ(
        host_target_node_freedom_signature(0), FreedomSignature::AllComponents
    );  // taget node has 6 DOFs which will be fixed
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_target_node_freedom_table(0, k), k + 17
        );  // target node DOFs: 17, 18, 19, 20, 21, 22
    }
}

TEST(TestCreateConstraintFreedomTable, SingleNodeConstraint_PrescribedBC) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto invalid_node = Node(0U, {0., 0., 0., 1., 0., 0., 0.});  // base node - index is 0
    auto node_1 = Node(1U, {1., 0., 0., 1., 0., 0., 0.});        // target node - index is 1
    auto prescribed_bc =
        constraints::Constraint(0, constraints::ConstraintType::PrescribedBC, {0, 1});
    auto constraints =
        Constraints<DeviceType>(std::array{prescribed_bc}, std::array{invalid_node, node_1});

    auto elements = Elements<DeviceType>();
    auto state = State<DeviceType>(2U);  // 2 nodes in the system
    assemble_node_freedom_allocation_table(state, elements, constraints);

    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 24UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    create_constraint_freedom_table(constraints, state);

    const auto host_target_node_index =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_index);
    const auto host_target_node_freedom_signature =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_freedom_signature);
    const auto host_target_node_freedom_table =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_freedom_table);

    EXPECT_EQ(host_target_node_index(0), 1);  // target node index is 4
    EXPECT_EQ(
        host_target_node_freedom_signature(0), FreedomSignature::AllComponents
    );  // taget node has 6 DOFs which will be prescribed
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_target_node_freedom_table(0, k), k + 24
        );  // target node DOFs: 24, 25, 26, 27, 28, 29
    }
}

TEST(TestCreateConstraintFreedomTable, DoubeNodeConstraint_RigidBC) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto node_1 = Node(0U, {0., 0., 0., 1., 0., 0., 0.});  // base node - index is 0
    auto node_2 = Node(1U, {1., 0., 0., 1., 0., 0., 0.});  // target node - index is 1
    auto rigid_bc = constraints::Constraint(0, constraints::ConstraintType::RigidJoint, {0, 1});
    auto constraints = Constraints<DeviceType>(std::array{rigid_bc}, std::array{node_1, node_2});

    auto elements = Elements<DeviceType>();
    auto state = State<DeviceType>(2U);  // 2 nodes in the system
    assemble_node_freedom_allocation_table(state, elements, constraints);

    constexpr auto host_node_freedom_map_table_data = std::array{0UL, 6UL};
    const auto host_node_freedom_map_table = Kokkos::View<size_t[2], Kokkos::HostSpace>::const_type(
        host_node_freedom_map_table_data.data()
    );
    const auto mirror_node_freedom_map_table =
        Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.node_freedom_map_table);
    Kokkos::deep_copy(mirror_node_freedom_map_table, host_node_freedom_map_table);
    Kokkos::deep_copy(state.node_freedom_map_table, mirror_node_freedom_map_table);

    create_constraint_freedom_table(constraints, state);

    const auto host_base_node_index =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.base_node_index);
    const auto host_target_node_index =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_index);
    const auto host_base_node_freedom_signature =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.base_node_freedom_signature);
    const auto host_target_node_freedom_signature =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_freedom_signature);
    const auto host_base_node_freedom_table =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.base_node_freedom_table);
    const auto host_target_node_freedom_table =
        create_mirror_view_and_copy(Kokkos::HostSpace(), constraints.target_node_freedom_table);

    EXPECT_EQ(host_base_node_index(0), 0);    // base node index is 0
    EXPECT_EQ(host_target_node_index(0), 1);  // target node index is 1
    EXPECT_EQ(
        host_base_node_freedom_signature(0), FreedomSignature::AllComponents
    );  // base node has 6 DOFs
    EXPECT_EQ(
        host_target_node_freedom_signature(0), FreedomSignature::AllComponents
    );  // target node has 6 DOFs which will be fixed relative to base node
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(host_base_node_freedom_table(0, k), k);  // base node DOFs: 0, 1, 2, 3, 4, 5
    }
    for (auto k : std::views::iota(0, 6)) {
        EXPECT_EQ(
            host_target_node_freedom_table(0, k), k + 6
        );  // target node DOFs: 6, 7, 8, 9, 10, 11
    }
}

}  // namespace kynema_fmb::dof::tests
