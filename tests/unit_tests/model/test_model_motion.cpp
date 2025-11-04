#include <gtest/gtest.h>

#include "model/model.hpp"

namespace kynema::tests {

class ModelMotionTest : public ::testing::Test {
protected:
    ModelMotionTest() {
        // node 1 is at origin
        node_1 = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
        // node 2 is at x=1
        node_2 = model.AddNode().SetPosition(1., 0., 0., 1., 0., 0., 0.).Build();
        // create a beam along the x-axis
        beam_element = model.AddBeamElement(std::array{node_1, node_2}, sections, quadrature);
    }

    Model model;
    size_t node_1;
    size_t node_2;
    size_t beam_element;
    std::vector<BeamSection> sections;
    std::vector<std::array<double, 2>> quadrature;
};

TEST_F(ModelMotionTest, TranslateBeam) {
    const auto displacement = std::array{0., 1., 2.};
    model.TranslateBeam(beam_element, displacement);

    // Check that nodes were translated
    auto n1 = model.GetNode(node_1);
    auto n2 = model.GetNode(node_2);

    // Node 1 should be translated by the displacement i.e. at (0, 1, 2)
    ASSERT_EQ(n1.x0[0], 0.);
    ASSERT_EQ(n1.x0[1], 1.);
    ASSERT_EQ(n1.x0[2], 2.);

    // Node 2 should be translated by the displacement i.e. at (1, 1, 2)
    ASSERT_EQ(n2.x0[0], 1.);
    ASSERT_EQ(n2.x0[1], 1.);
    ASSERT_EQ(n2.x0[2], 2.);

    // Test additional translation (should accumulate)
    const auto displacement2 = std::array{1., 1., -1.};
    model.TranslateBeam(beam_element, displacement2);

    // Check updated positions
    n1 = model.GetNode(node_1);
    n2 = model.GetNode(node_2);

    // Node 1 should be translated by the displacement i.e. at (1, 2, 1)
    ASSERT_EQ(n1.x0[0], 1.);
    ASSERT_EQ(n1.x0[1], 2.);
    ASSERT_EQ(n1.x0[2], 1.);

    // Node 2 should be translated by the displacement i.e. at (2, 2, 1)
    ASSERT_EQ(n2.x0[0], 2.);
    ASSERT_EQ(n2.x0[1], 2.);
    ASSERT_EQ(n2.x0[2], 1.);
}

TEST_F(ModelMotionTest, RotateBeamAboutOrigin) {
    // Rotate 90 degrees around Z axis using quaternion [cos(π/4), 0, 0, sin(π/4)]
    const auto rotation = std::array{0.7071068, 0., 0., 0.7071068};  // cos(45°), 0, 0, sin(45°)
    const auto center = std::array{0., 0., 0.};
    model.RotateBeamAboutPoint(beam_element, rotation, center);

    // Check that nodes were rotated properly
    auto n1 = model.GetNode(node_1);
    auto n2 = model.GetNode(node_2);

    // Node 1 is at origin, so it shouldn't move
    ASSERT_NEAR(n1.x0[0], 0., 1e-6);
    ASSERT_NEAR(n1.x0[1], 0., 1e-6);
    ASSERT_NEAR(n1.x0[2], 0., 1e-6);

    // Node 2 should rotate 90 degrees around the z-axis
    // Final position should be approximately (0, 1, 0)
    ASSERT_NEAR(n2.x0[0], 0., 1e-6);
    ASSERT_NEAR(n2.x0[1], 1., 1e-6);
    ASSERT_NEAR(n2.x0[2], 0., 1e-6);
}

TEST_F(ModelMotionTest, RotateBeamAboutArbitraryPoint) {
    // First translate the beam to have Node 1 at (1, 1, 0) and Node 2 at (2, 1, 0)
    model.TranslateBeam(beam_element, std::array{1., 1., 0.});

    // Now rotate around point (1, 1, 0) by 90 degrees around y-axis
    const auto rotation = std::array{0.7071068, 0., 0.7071068, 0.};  // 90 degrees around y-axis
    const auto center = std::array{1., 1., 0.};
    model.RotateBeamAboutPoint(beam_element, rotation, center);

    // Check that nodes were rotated properly around the new center
    auto n1 = model.GetNode(node_1);
    auto n2 = model.GetNode(node_2);

    // Node 1 was at the center of rotation, so it shouldn't move
    ASSERT_NEAR(n1.x0[0], 1., 1e-6);
    ASSERT_NEAR(n1.x0[1], 1., 1e-6);
    ASSERT_NEAR(n1.x0[2], 0., 1e-6);

    // Node 2 was at (2, 1, 0) and should rotate around y to (1, 1, -1)
    ASSERT_NEAR(n2.x0[0], 1., 1e-6);
    ASSERT_NEAR(n2.x0[1], 1., 1e-6);
    ASSERT_NEAR(n2.x0[2], -1., 1e-6);
}

TEST_F(ModelMotionTest, SetBeamVelocityAboutOrigin) {
    // Set a velocity about the origin
    const auto velocity = std::array{
        1., 0., 0.,  // Linear velocity in x
        0., 0., 1.   // Angular velocity around z
    };
    const auto point = std::array{0., 0., 0.};
    model.SetBeamVelocityAboutPoint(beam_element, velocity, point);

    // Check velocity at first node (at origin)
    auto n1 = model.GetNode(node_1);
    ASSERT_EQ(n1.v[0], 1.);  // Linear velocity in x
    ASSERT_EQ(n1.v[1], 0.);
    ASSERT_EQ(n1.v[2], 0.);
    ASSERT_EQ(n1.v[3], 0.);  // Angular velocity
    ASSERT_EQ(n1.v[4], 0.);
    ASSERT_EQ(n1.v[5], 1.);

    // Check velocity at second node (at x=1)
    // For angular velocity around z, this should add a y-component to velocity
    auto n2 = model.GetNode(node_2);
    ASSERT_EQ(n2.v[0], 1.);  // Linear velocity in x
    ASSERT_EQ(n2.v[1], 1.);  // Additional y velocity from rotation
    ASSERT_EQ(n2.v[2], 0.);
    ASSERT_EQ(n2.v[3], 0.);  // Angular velocity
    ASSERT_EQ(n2.v[4], 0.);
    ASSERT_EQ(n2.v[5], 1.);
}

TEST_F(ModelMotionTest, SetBeamVelocityAboutArbitraryPoint) {
    // Test with reference point at (1, 0, 0)
    const auto velocity = std::array{
        0., 1., 0.,  // Linear velocity in y
        0., 0., 1.   // Angular velocity around z
    };
    const auto point = std::array{1., 0., 0.};  // Reference point at node 2
    model.SetBeamVelocityAboutPoint(beam_element, velocity, point);

    // Check velocities
    auto n1 = model.GetNode(node_1);
    auto n2 = model.GetNode(node_2);

    // Node 1 is at (-1,0,0) relative to the reference point (1, 0, 0)
    // So it should get additional -1 in x direction from rotation around z
    ASSERT_EQ(n1.v[0], 0.);
    ASSERT_EQ(n1.v[1], 1. - 1.);  // Base y velocity + contribution from rotation = 0
    ASSERT_EQ(n1.v[2], 0.);
    ASSERT_EQ(n1.v[3], 0.);
    ASSERT_EQ(n1.v[4], 0.);
    ASSERT_EQ(n1.v[5], 1.);

    // Node 2 is at the reference point, so no additional velocity
    ASSERT_EQ(n2.v[0], 0.);
    ASSERT_EQ(n2.v[1], 1.);
    ASSERT_EQ(n2.v[2], 0.);
    ASSERT_EQ(n2.v[3], 0.);
    ASSERT_EQ(n2.v[4], 0.);
    ASSERT_EQ(n2.v[5], 1.);
}

TEST_F(ModelMotionTest, SetBeamAccelerationAboutOrigin) {
    // Set an acceleration about the origin
    const auto acceleration = std::array{
        1., 0., 0.,  // Linear acceleration in x
        0., 0., 1.   // Angular acceleration around z
    };
    const auto omega = std::array{0., 0., 1.};  // Angular velocity vector (1 rad/s around z-axis)
    const auto point = std::array{0., 0., 0.};  // Reference point at origin
    model.SetBeamAccelerationAboutPoint(beam_element, acceleration, omega, point);

    // Check acceleration at first node (at origin)
    auto n1 = model.GetNode(node_1);
    ASSERT_EQ(n1.vd[0], 1.);  // Linear acceleration in x
    ASSERT_EQ(n1.vd[1], 0.);
    ASSERT_EQ(n1.vd[2], 0.);
    ASSERT_EQ(n1.vd[3], 0.);
    ASSERT_EQ(n1.vd[4], 0.);
    ASSERT_EQ(n1.vd[5], 1.);  // Angular acceleration around z

    // Check acceleration at second node (at x=1)
    // For node at (1,0,0):
    // - Angular acceleration contribution: α × r = (0,0,1) × (1,0,0) = (0,1,0)
    // - Centripetal acceleration: ω × (ω × r) = (0,0,1) × ((0,0,1) × (1,0,0))
    //                                         = (0,0,1) × (0,1,0) = (-1,0,0)
    auto n2 = model.GetNode(node_2);
    ASSERT_EQ(n2.vd[0], 1. - 1.);  // Linear + centripetal acceleration (1 - 1 = 0)
    ASSERT_EQ(n2.vd[1], 1.);       // Angular acceleration contribution
    ASSERT_EQ(n2.vd[2], 0.);
    ASSERT_EQ(n2.vd[3], 0.);  // Angular acceleration
    ASSERT_EQ(n2.vd[4], 0.);
    ASSERT_EQ(n2.vd[5], 1.);
}

TEST_F(ModelMotionTest, SetBeamAccelerationAboutArbitraryPoint) {
    // Test with reference point at node 2
    const auto acceleration = std::array{
        0., 1., 0.,  // Linear acceleration in y
        0., 1., 0.   // Angular acceleration around y
    };
    const auto omega = std::array{0., 1., 0.};  // Angular velocity vector (1 rad/s around y-axis)
    const auto point = std::array{1., 0., 0.};  // Reference point at node 2
    model.SetBeamAccelerationAboutPoint(beam_element, acceleration, omega, point);

    // Check updated accelerations
    auto n1 = model.GetNode(node_1);
    auto n2 = model.GetNode(node_2);

    // Node 1 is at position (-1,0,0) relative to reference point (1,0,0)
    // - Angular acceleration contribution: α × r = (0,1,0) × (-1,0,0) = (0,0,1)
    // - Centripetal acceleration: ω × (ω × r) = (0,1,0) × ((0,1,0) × (-1,0,0))
    //                                         = (0,1,0) × (0,0,1) = (1,0,0)
    ASSERT_EQ(n1.vd[0], 0. + 1.);  // Linear + centripetal
    ASSERT_EQ(n1.vd[1], 1.);
    ASSERT_EQ(n1.vd[2], 1.);  // From angular acceleration
    ASSERT_EQ(n1.vd[3], 0.);
    ASSERT_EQ(n1.vd[4], 1.);
    ASSERT_EQ(n1.vd[5], 0.);

    // Node2 is at the reference point, so no additional acceleration from rotation
    ASSERT_EQ(n2.vd[0], 0.);
    ASSERT_EQ(n2.vd[1], 1.);
    ASSERT_EQ(n2.vd[2], 0.);
    ASSERT_EQ(n2.vd[3], 0.);
    ASSERT_EQ(n2.vd[4], 1.);
    ASSERT_EQ(n2.vd[5], 0.);
}

TEST_F(ModelMotionTest, ComplexMotionSequence) {
    // This test applies a sequence of motions to verify correct composition

    // ----------------------------------------------
    // Step 1: Translate beam up along z-axis
    // ----------------------------------------------
    model.TranslateBeam(beam_element, std::array{0., 0., 1.});

    // ----------------------------------------------
    // Step 2: Rotate 90 degrees around x-axis
    // ----------------------------------------------
    const auto rotation_x = std::array{0.7071068, 0.7071068, 0., 0.};  // 90 degrees around x-axis
    const auto origin = std::array{0., 0., 0.};
    model.RotateBeamAboutPoint(beam_element, rotation_x, origin);

    // Check intermediate positions after translation and rotation
    auto n1 = model.GetNode(node_1);
    auto n2 = model.GetNode(node_2);

    // Node 1 should be at (0,0,1) rotated to (0,-1,0)
    ASSERT_NEAR(n1.x0[0], 0., 1e-6);
    ASSERT_NEAR(n1.x0[1], -1., 1e-6);
    ASSERT_NEAR(n1.x0[2], 0., 1e-6);

    // Node 2 should be at (1,0,1) rotated to (1,-1,0)
    ASSERT_NEAR(n2.x0[0], 1., 1e-6);
    ASSERT_NEAR(n2.x0[1], -1., 1e-6);
    ASSERT_NEAR(n2.x0[2], 0., 1e-6);

    // ----------------------------------------------
    // Step 3: Set velocity about the origin
    // ----------------------------------------------
    // Linear velocity in y direction, angular velocity around z
    const auto velocity = std::array{0., 2., 0., 0., 0., 1.};
    model.SetBeamVelocityAboutPoint(beam_element, velocity, origin);

    // Check velocities after applying step 3
    n1 = model.GetNode(node_1);
    n2 = model.GetNode(node_2);

    // Node 1 velocity: linear + angular contribution
    ASSERT_NEAR(n1.v[0], 1., 1e-6);  // Angular contribution
    ASSERT_NEAR(n1.v[1], 2., 1e-6);  // Linear y velocity
    ASSERT_NEAR(n1.v[2], 0., 1e-6);
    ASSERT_NEAR(n1.v[3], 0., 1e-6);
    ASSERT_NEAR(n1.v[4], 0., 1e-6);
    ASSERT_NEAR(n1.v[5], 1., 1e-6);  // Angular z velocity

    // Node 2 velocity: should have additional component from x-position
    ASSERT_NEAR(n2.v[0], 1., 1e-6);  // Angular contribution
    ASSERT_NEAR(n2.v[1], 3., 1e-6);  // Linear y velocity + Angular contribution
    ASSERT_NEAR(n2.v[2], 0., 1e-6);
    ASSERT_NEAR(n2.v[3], 0., 1e-6);
    ASSERT_NEAR(n2.v[4], 0., 1e-6);
    ASSERT_NEAR(n2.v[5], 1., 1e-6);  // Angular z velocity

    // ----------------------------------------------
    // Step 4: Set acceleration about a different point (the first node)
    // ----------------------------------------------
    const auto acceleration = std::array{
        1., 0., 0.,  // Linear x accleration
        0., 1., 0.   // Angular y accleration
    };
    const auto omega = std::array{0., 0., 1.};       // Current angular velocity is around z
    const auto ref_point = std::array{0., -1., 0.};  // Position of first node after transformations
    model.SetBeamAccelerationAboutPoint(beam_element, acceleration, omega, ref_point);

    // Check final accelerations
    n1 = model.GetNode(node_1);
    n2 = model.GetNode(node_2);

    // Node 1 is at reference point, so it just gets the linear and angular components
    ASSERT_NEAR(n1.vd[0], 1., 1e-6);  // Linear x acceleration
    ASSERT_NEAR(n1.vd[1], 0., 1e-6);  // Angular y acceleration
    ASSERT_NEAR(n1.vd[2], 0., 1e-6);
    ASSERT_NEAR(n1.vd[3], 0., 1e-6);
    ASSERT_NEAR(n1.vd[4], 1., 1e-6);  // Angular y acceleration
    ASSERT_NEAR(n1.vd[5], 0., 1e-6);

    // Node 2 is at (1, -1, 0) relative to (0, -1, 0), so (1, 0, 0) displacement vector
    // - Angular acceln contribution: (0,1,0) × (1,0,0) = (0,0,-1)
    // - Centripetal acceln from current omega: (0,0,1) × ((0,0,1) × (1,0,0))
    //                                          = (0,0,1) × (0,1,0) = (-1,0,0)
    ASSERT_NEAR(n2.vd[0], 1. - 1., 1e-6);  // Linear + centripetal = 0
    ASSERT_NEAR(n2.vd[1], 0., 1e-6);       // Angular y acceleration
    ASSERT_NEAR(n2.vd[2], -1., 1e-6);      // From angular acceleration contribution
    ASSERT_NEAR(n2.vd[3], 0., 1e-6);
    ASSERT_NEAR(n2.vd[4], 1., 1e-6);  // Angular y acceleration
    ASSERT_NEAR(n2.vd[5], 0., 1e-6);
}

}  // namespace kynema::tests
