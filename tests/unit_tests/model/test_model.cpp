#include <gtest/gtest.h>

#include "model/model.hpp"

namespace kynema::tests {

TEST(ModelTest, AddNodeToModel) {
    Model model;
    ASSERT_EQ(model.NumNodes(), 0);

    // Add a node to the model
    constexpr auto pos = std::array{0., 0., 0.};
    constexpr auto rot = std::array{1., 0., 0., 0.};
    constexpr auto v = std::array{0., 0., 0.};
    constexpr auto omega = std::array{0., 0., 0.};
    auto node_id = model.AddNode()
                       .SetPosition(pos[0], pos[1], pos[2], rot[0], rot[1], rot[2], rot[3])
                       .SetVelocity(v[0], v[1], v[2], omega[0], omega[1], omega[2])
                       .Build();
    auto node = model.GetNode(node_id);
    ASSERT_EQ(node_id, 0);
    ASSERT_EQ(node.id, 0);
    ASSERT_EQ(model.NumNodes(), 1);

    auto nodes = model.GetNodes();
    ASSERT_EQ(nodes.size(), 1);
}

TEST(ModelTest, AddBeamElementToModel) {
    Model model;

    constexpr auto pos = std::array{0., 0., 0.};
    constexpr auto rot = std::array{1., 0., 0., 0.};
    constexpr auto v = std::array{0., 0., 0.};
    constexpr auto omega = std::array{0., 0., 0.};

    // Add couple of nodes to the model
    auto node1_id = model.AddNode()
                        .SetElemLocation(0.)
                        .SetPosition(pos[0], pos[1], pos[2], rot[0], rot[1], rot[2], rot[3])
                        .SetVelocity(v[0], v[1], v[2], omega[0], omega[1], omega[2])
                        .Build();

    auto node2_id = model.AddNode()
                        .SetElemLocation(1.)
                        .SetPosition(pos[0], pos[1], pos[2], rot[0], rot[1], rot[2], rot[3])
                        .SetVelocity(v[0], v[1], v[2], omega[0], omega[1], omega[2])
                        .Build();

    // Add a beam element to the model
    constexpr auto sections = std::array<BeamSection, 0>{};
    constexpr auto quadrature = std::array<std::array<double, 2>, 0>{};
    model.AddBeamElement(std::array{node1_id, node2_id}, sections, quadrature);

    ASSERT_EQ(model.NumBeamElements(), 1);

    auto elements = model.GetBeamElements();
    ASSERT_EQ(elements.size(), 1);
}

TEST(ModelTest, AddMassElementToModel) {
    Model model;
    auto node_id = model.AddNode().Build();
    auto mass_matrix = std::array<std::array<double, 6>, 6>{};
    model.AddMassElement(node_id, mass_matrix);

    ASSERT_EQ(model.NumMassElements(), 1);

    auto elements = model.GetMassElements();
    ASSERT_EQ(elements.size(), 1);
}

TEST(ModelTest, ModelConstructorWithDefaults) {
    const Model model;
    ASSERT_EQ(model.NumNodes(), 0);
    ASSERT_EQ(model.NumBeamElements(), 0);
    ASSERT_EQ(model.NumMassElements(), 0);
    ASSERT_EQ(model.NumConstraints(), 0);
}

TEST(ModelTest, ModelCreateState) {
    Model model;

    // Rotation of 1 radian around x
    const auto R1 =
        Eigen::Quaternion<double>(Eigen::AngleAxis<double>(1., Eigen::Matrix<double, 3, 1>::Unit(0))
        );
    const auto R2 =
        Eigen::Quaternion<double>(Eigen::AngleAxis<double>(1., Eigen::Matrix<double, 3, 1>::Unit(1))
        );

    // Create node with initial position and displacement from initial position
    static_cast<void>(model.AddNode()
                          .SetPosition(1., 2., 3., R1.w(), R1.x(), R1.y(), R1.z())
                          .SetDisplacement(3., 2., 1., R2.w(), R2.x(), R2.y(), R2.z())
                          .Build());

    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    // Create state object from model
    auto state = model.CreateState<DeviceType>();

    // Verify initial position
    const auto x0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.x0);
    const auto exact_x0 = std::array{1., 2., 3., R1.w(), R1.x(), R1.y(), R1.z()};
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(x0(0, i), exact_x0[i], 1.e-15);
    }

    // Verify initial displacement
    const auto q = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.q);
    const auto exact_q = std::array{3., 2., 1., R2.w(), R2.x(), R2.y(), R2.z()};
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(q(0, i), exact_q[i], 1.e-15);
    }

    // Verify current position (initial position plus displacement)
    const auto Rt = R2 * R1;
    const auto x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.x);
    const auto exact_x = std::array{4., 4., 4., Rt.w(), Rt.x(), Rt.y(), Rt.z()};
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(x(0, i), exact_x[i], 1.e-15);
    }
}

TEST(ModelTest, ModelCreateSystem) {
    Model model;

    // Rotation of 1 radian around x
    const auto R1 =
        Eigen::Quaternion<double>(Eigen::AngleAxis<double>(1., Eigen::Matrix<double, 3, 1>::Unit(0))
        );
    const auto R2 =
        Eigen::Quaternion<double>(Eigen::AngleAxis<double>(1., Eigen::Matrix<double, 3, 1>::Unit(1))
        );

    // Create node with initial position and displacement from initial position
    static_cast<void>(model.AddNode()
                          .SetPosition(1., 2., 3., R1.w(), R1.x(), R1.y(), R1.z())
                          .SetDisplacement(3., 2., 1., R2.w(), R2.x(), R2.y(), R2.z())
                          .Build());

    // Create state object from model
    auto [state, elements, constraints] = model.CreateSystem();

    // Verify initial position
    const auto x0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.x0);
    const auto exact_x0 = std::array{1., 2., 3., R1.w(), R1.x(), R1.y(), R1.z()};
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(x0(0, i), exact_x0[i], 1.e-15);
    }

    // Verify initial displacement
    const auto q = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.q);
    const auto exact_q = std::array{3., 2., 1., R2.w(), R2.x(), R2.y(), R2.z()};
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(q(0, i), exact_q[i], 1.e-15);
    }

    // Verify current position (initial position plus displacement)
    const auto Rt = R2 * R1;
    const auto x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.x);
    const auto exact_x = std::array{4., 4., 4., Rt.w(), Rt.x(), Rt.y(), Rt.z()};
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(x(0, i), exact_x[i], 1.e-15);
    }

    EXPECT_EQ(elements.NumElementsInSystem(), 0);
    EXPECT_EQ(constraints.num_constraints, 0);
}

}  // namespace kynema::tests
