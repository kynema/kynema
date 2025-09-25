#include <array>
#include <cmath>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "interfaces/components/beam.hpp"
#include "interfaces/components/beam_input.hpp"
#include "model/model.hpp"

namespace kynema::tests {

TEST(BeamComponentTest, InitialBeamHasCorrectRotation) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.coordinates =
        std::vector{std::array{0., 0., 0.}, std::array{0.5, 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(0.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto& beam_node_ids = model.GetBeamElement(0).node_ids;

    const auto& root_node = model.GetNode(beam_node_ids[0]);
    const auto root_rotation_quaternion =
        std::array{root_node.x0[3], root_node.x0[4], root_node.x0[5], root_node.x0[6]};

    const auto root_rotation_matrix = math::QuaternionToRotationMatrix(root_rotation_quaternion);

    EXPECT_NEAR(root_node.x0[0], 0., 1.e-16);
    EXPECT_NEAR(root_node.x0[1], 0., 1.e-16);
    EXPECT_NEAR(root_node.x0[2], 0., 1.e-16);

    EXPECT_NEAR(root_rotation_matrix[0][0], 1., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[0][1], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[0][2], 0., 1.e-16);

    EXPECT_NEAR(root_rotation_matrix[1][0], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[1][1], 1., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[1][2], 0., 1.e-16);

    EXPECT_NEAR(root_rotation_matrix[2][0], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[2][1], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[2][2], 1., 1.e-16);

    const auto& middle_node = model.GetNode(beam_node_ids[1]);
    const auto middle_rotation_quaternion =
        std::array{middle_node.x0[3], middle_node.x0[4], middle_node.x0[5], middle_node.x0[6]};

    const auto middle_rotation_matrix = math::QuaternionToRotationMatrix(middle_rotation_quaternion);

    EXPECT_NEAR(middle_node.x0[0], .5, 1.e-16);
    EXPECT_NEAR(middle_node.x0[1], 0., 1.e-16);
    EXPECT_NEAR(middle_node.x0[2], 0., 1.e-16);

    EXPECT_NEAR(middle_rotation_matrix[0][0], 1., 1.e-16);
    EXPECT_NEAR(middle_rotation_matrix[0][1], 0., 1.e-16);
    EXPECT_NEAR(middle_rotation_matrix[0][2], 0., 1.e-16);

    EXPECT_NEAR(middle_rotation_matrix[1][0], 0., 1.e-16);
    EXPECT_NEAR(middle_rotation_matrix[1][1], 1., 1.e-16);
    EXPECT_NEAR(middle_rotation_matrix[1][2], 0., 1.e-16);

    EXPECT_NEAR(middle_rotation_matrix[2][0], 0., 1.e-16);
    EXPECT_NEAR(middle_rotation_matrix[2][1], 0., 1.e-16);
    EXPECT_NEAR(middle_rotation_matrix[2][2], 1., 1.e-16);

    const auto& end_node = model.GetNode(beam_node_ids[2]);
    const auto end_rotation_quaternion =
        std::array{end_node.x0[3], end_node.x0[4], end_node.x0[5], end_node.x0[6]};

    const auto end_rotation_matrix = math::QuaternionToRotationMatrix(end_rotation_quaternion);

    EXPECT_NEAR(end_node.x0[0], 1., 1.e-16);
    EXPECT_NEAR(end_node.x0[1], 0., 1.e-16);
    EXPECT_NEAR(end_node.x0[2], 0., 1.e-16);

    EXPECT_NEAR(end_rotation_matrix[0][0], 1., 1.e-16);
    EXPECT_NEAR(end_rotation_matrix[0][1], 0., 1.e-16);
    EXPECT_NEAR(end_rotation_matrix[0][2], 0., 1.e-16);

    EXPECT_NEAR(end_rotation_matrix[1][0], 0., 1.e-16);
    EXPECT_NEAR(end_rotation_matrix[1][1], 1., 1.e-16);
    EXPECT_NEAR(end_rotation_matrix[1][2], 0., 1.e-16);

    EXPECT_NEAR(end_rotation_matrix[2][0], 0., 1.e-16);
    EXPECT_NEAR(end_rotation_matrix[2][1], 0., 1.e-16);
    EXPECT_NEAR(end_rotation_matrix[2][2], 1., 1.e-16);
}

TEST(BeamComponentTest, UnrotatedBeamHasIdentityRotationMatrix) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.coordinates =
        std::vector{std::array{0., 0., 0.}, std::array{0.5, 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(0.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto& beam_node_ids = model.GetBeamElement(0).node_ids;

    const auto& root_node = model.GetNode(beam_node_ids[0]);
    const auto root_rotation_quaternion =
        std::array{root_node.x0[3], root_node.x0[4], root_node.x0[5], root_node.x0[6]};

    const auto root_rotation_matrix = math::QuaternionToRotationMatrix(root_rotation_quaternion);

    EXPECT_NEAR(root_rotation_matrix[0][0], 1., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[0][1], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[0][2], 0., 1.e-16);

    EXPECT_NEAR(root_rotation_matrix[1][0], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[1][1], 1., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[1][2], 0., 1.e-16);

    EXPECT_NEAR(root_rotation_matrix[2][0], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[2][1], 0., 1.e-16);
    EXPECT_NEAR(root_rotation_matrix[2][2], 1., 1.e-16);
}

TEST(BeamComponentTest, RotatedBeamAboutYAxisPointsAlongZAxis) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.coordinates =
        std::vector{std::array{0., 0., 0.}, std::array{0.5, 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(0.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.root.position =
        std::array{0., 0., 0., std::cos(std::numbers::pi / 4.), 0., std::sin(std::numbers::pi / 4.),
                   0.};

    auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto& beam_node_ids = model.GetBeamElement(0).node_ids;

    const auto& root_node = model.GetNode(beam_node_ids[0]);
    const auto root_rotation_quaternion =
        std::array{root_node.x0[3], root_node.x0[4], root_node.x0[5], root_node.x0[6]};

    const auto root_rotation_matrix = math::QuaternionToRotationMatrix(root_rotation_quaternion);

    EXPECT_NEAR(root_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(root_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(root_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[0][0], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[0][1], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[0][2], 1., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[1][0], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[1][1], 1., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[1][2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[2][0], -1., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[2][1], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[2][2], 0., 1.e-15);

    const auto& middle_node = model.GetNode(beam_node_ids[1]);
    const auto middle_rotation_quaternion =
        std::array{middle_node.x0[3], middle_node.x0[4], middle_node.x0[5], middle_node.x0[6]};

    const auto middle_rotation_matrix = math::QuaternionToRotationMatrix(middle_rotation_quaternion);

    EXPECT_NEAR(middle_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(middle_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(middle_node.x0[2], -.5, 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[0][0], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[0][1], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[0][2], 1., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[1][0], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[1][1], 1., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[1][2], 0., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[2][0], -1., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[2][1], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[2][2], 0., 1.e-15);

    const auto& end_node = model.GetNode(beam_node_ids[2]);
    const auto end_rotation_quaternion =
        std::array{end_node.x0[3], end_node.x0[4], end_node.x0[5], end_node.x0[6]};

    const auto end_rotation_matrix = math::QuaternionToRotationMatrix(end_rotation_quaternion);

    EXPECT_NEAR(end_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(end_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(end_node.x0[2], -1., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[0][0], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[0][1], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[0][2], 1., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[1][0], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[1][1], 1., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[1][2], 0., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[2][0], -1., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[2][1], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[2][2], 0., 1.e-15);
}

TEST(BeamComponentTest, RotatedBeamAboutZAxisPointsAlongYAxis) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.coordinates =
        std::vector{std::array{0., 0., 0.}, std::array{0.5, 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(0.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.root.position = std::array{
        0., 0., 0., std::cos(std::numbers::pi / 4.), 0., 0., std::sin(std::numbers::pi / 4.)
    };

    auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto& beam_node_ids = model.GetBeamElement(0).node_ids;

    const auto& root_node = model.GetNode(beam_node_ids[0]);
    const auto root_rotation_quaternion =
        std::array{root_node.x0[3], root_node.x0[4], root_node.x0[5], root_node.x0[6]};

    const auto root_rotation_matrix = math::QuaternionToRotationMatrix(root_rotation_quaternion);

    EXPECT_NEAR(root_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(root_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(root_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[0][0], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[0][1], -1., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[0][2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[1][0], 1., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[1][1], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[1][2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[2][0], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[2][1], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[2][2], 1., 1.e-15);

    const auto& middle_node = model.GetNode(beam_node_ids[1]);
    const auto middle_rotation_quaternion =
        std::array{middle_node.x0[3], middle_node.x0[4], middle_node.x0[5], middle_node.x0[6]};

    const auto middle_rotation_matrix = math::QuaternionToRotationMatrix(middle_rotation_quaternion);

    EXPECT_NEAR(middle_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(middle_node.x0[1], .5, 1.e-15);
    EXPECT_NEAR(middle_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[0][0], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[0][1], -1., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[0][2], 0., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[1][0], 1., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[1][1], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[1][2], 0., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[2][0], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[2][1], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[2][2], 1., 1.e-15);

    const auto& end_node = model.GetNode(beam_node_ids[2]);
    const auto end_rotation_quaternion =
        std::array{end_node.x0[3], end_node.x0[4], end_node.x0[5], end_node.x0[6]};

    const auto end_rotation_matrix = math::QuaternionToRotationMatrix(end_rotation_quaternion);

    EXPECT_NEAR(end_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(end_node.x0[1], 1., 1.e-15);
    EXPECT_NEAR(end_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[0][0], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[0][1], -1., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[0][2], 0., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[1][0], 1., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[1][1], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[1][2], 0., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[2][0], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[2][1], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[2][2], 1., 1.e-15);
}

TEST(BeamComponentTest, RotatedBeamAboutXAxisStillPointsAlongXAxis) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.coordinates =
        std::vector{std::array{0., 0., 0.}, std::array{0.5, 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 0.5, 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(0.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.root.position =
        std::array{0., 0., 0., std::cos(std::numbers::pi / 4.), std::sin(std::numbers::pi / 4.),
                   0., 0.};

    auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto& beam_node_ids = model.GetBeamElement(0).node_ids;

    const auto& root_node = model.GetNode(beam_node_ids[0]);
    const auto root_rotation_quaternion =
        std::array{root_node.x0[3], root_node.x0[4], root_node.x0[5], root_node.x0[6]};

    const auto root_rotation_matrix = math::QuaternionToRotationMatrix(root_rotation_quaternion);

    EXPECT_NEAR(root_node.x0[0], 0., 1.e-15);
    EXPECT_NEAR(root_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(root_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[0][0], 1., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[0][1], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[0][2], 0., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[1][0], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[1][1], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[1][2], -1., 1.e-15);

    EXPECT_NEAR(root_rotation_matrix[2][0], 0., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[2][1], 1., 1.e-15);
    EXPECT_NEAR(root_rotation_matrix[2][2], 0., 1.e-15);

    const auto& middle_node = model.GetNode(beam_node_ids[1]);
    const auto middle_rotation_quaternion =
        std::array{middle_node.x0[3], middle_node.x0[4], middle_node.x0[5], middle_node.x0[6]};

    const auto middle_rotation_matrix = math::QuaternionToRotationMatrix(middle_rotation_quaternion);

    EXPECT_NEAR(middle_node.x0[0], .5, 1.e-15);
    EXPECT_NEAR(middle_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(middle_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[0][0], 1., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[0][1], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[0][2], 0., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[1][0], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[1][1], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[1][2], -1., 1.e-15);

    EXPECT_NEAR(middle_rotation_matrix[2][0], 0., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[2][1], 1., 1.e-15);
    EXPECT_NEAR(middle_rotation_matrix[2][2], 0., 1.e-15);

    const auto& end_node = model.GetNode(beam_node_ids[2]);
    const auto end_rotation_quaternion =
        std::array{end_node.x0[3], end_node.x0[4], end_node.x0[5], end_node.x0[6]};

    const auto end_rotation_matrix = math::QuaternionToRotationMatrix(end_rotation_quaternion);

    EXPECT_NEAR(end_node.x0[0], 1., 1.e-15);
    EXPECT_NEAR(end_node.x0[1], 0., 1.e-15);
    EXPECT_NEAR(end_node.x0[2], 0., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[0][0], 1., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[0][1], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[0][2], 0., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[1][0], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[1][1], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[1][2], -1., 1.e-15);

    EXPECT_NEAR(end_rotation_matrix[2][0], 0., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[2][1], 1., 1.e-15);
    EXPECT_NEAR(end_rotation_matrix[2][2], 0., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_TwoSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 1., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_ThreeSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[1][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 1., 1.e-15);

    EXPECT_NEAR(quadrature[2][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .5, 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_FourSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.75, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[1][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], .75, 1.e-15);

    EXPECT_NEAR(quadrature[2][0], .5, 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[3][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[3][1], .25, 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_FourSections_NoRefinement_Length2) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{2., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.75, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[1][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], .75, 1.e-15);

    EXPECT_NEAR(quadrature[2][0], .5, 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[3][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[3][1], .25, 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_ElevenSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.1, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.2, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.3, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.4, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.6, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.7, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.8, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.9, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -.8, 1.e-15);
    EXPECT_NEAR(quadrature[1][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[2][0], -.6, 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[3][0], -.4, 1.e-15);
    EXPECT_NEAR(quadrature[3][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[4][0], -.2, 1.e-15);
    EXPECT_NEAR(quadrature[4][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[5][0], -0., 1.e-15);
    EXPECT_NEAR(quadrature[5][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[6][0], .2, 1.e-15);
    EXPECT_NEAR(quadrature[6][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[7][0], .4, 1.e-15);
    EXPECT_NEAR(quadrature[7][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[8][0], .6, 1.e-15);
    EXPECT_NEAR(quadrature[8][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[9][0], .8, 1.e-15);
    EXPECT_NEAR(quadrature[9][1], .2, 1.e-15);

    EXPECT_NEAR(quadrature[10][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[10][1], .1, 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_TwoSections_OneRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1. / 3., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 4. / 3., 1.e-15);

    EXPECT_NEAR(quadrature[2][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[2][1], 1. / 3., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_TwoSections_TwoRefinements) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 2UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1. / 6., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -std::sqrt(1. / 5.), 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 5. / 6., 1.e-15);

    EXPECT_NEAR(quadrature[2][0], std::sqrt(1. / 5.), 1.e-15);
    EXPECT_NEAR(quadrature[2][1], 5. / 6., 1.e-15);

    EXPECT_NEAR(quadrature[3][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[3][1], 1. / 6., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_ThreeSections_OneRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1. / 6., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -.5, 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 4. / 6., 1.e-15);

    EXPECT_NEAR(quadrature[2][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[2][1], 1. / 3., 1.e-15);

    EXPECT_NEAR(quadrature[3][0], .5, 1.e-15);
    EXPECT_NEAR(quadrature[3][1], 4. / 6., 1.e-15);

    EXPECT_NEAR(quadrature[4][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[4][1], 1. / 6., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_FourSections_TwoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.75, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 2UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1. / 12., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -.5 - (std::sqrt(1. / 5.) / 2.), 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 5. / 12., 1.e-15);

    EXPECT_NEAR(quadrature[2][0], -.5 + (std::sqrt(1. / 5.) / 2.), 1.e-15);
    EXPECT_NEAR(quadrature[2][1], 5. / 12., 1.e-15);

    EXPECT_NEAR(quadrature[3][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[3][1], (1. / 12.) + (1. / 24.), 1.e-15);

    EXPECT_NEAR(quadrature[4][0], .25 - (std::sqrt(1. / 5.) / 4.), 1.e-15);
    EXPECT_NEAR(quadrature[4][1], 5. / 24., 1.e-15);

    EXPECT_NEAR(quadrature[5][0], .25 + (std::sqrt(1. / 5.) / 4.), 1.e-15);
    EXPECT_NEAR(quadrature[5][1], 5. / 24., 1.e-15);

    EXPECT_NEAR(quadrature[6][0], .5, 1.e-15);
    EXPECT_NEAR(quadrature[6][1], 1. / 12., 1.e-15);

    EXPECT_NEAR(quadrature[7][0], .75 - (std::sqrt(1. / 5.) / 4.), 1.e-15);
    EXPECT_NEAR(quadrature[7][1], 5. / 24., 1.e-15);

    EXPECT_NEAR(quadrature[8][0], .75 + (std::sqrt(1. / 5.) / 4.), 1.e-15);
    EXPECT_NEAR(quadrature[8][1], 5. / 24., 1.e-15);

    EXPECT_NEAR(quadrature[9][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[9][1], 1. / 24., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GLL_ElevenSections_OneRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.1, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.2, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.3, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.4, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.6, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.7, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.8, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.9, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], (1. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -.9, 1.e-15);
    EXPECT_NEAR(quadrature[1][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[2][0], -.8, 1.e-15);
    EXPECT_NEAR(quadrature[2][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[3][0], -.7, 1.e-15);
    EXPECT_NEAR(quadrature[3][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[4][0], -.6, 1.e-15);
    EXPECT_NEAR(quadrature[4][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[5][0], -.5, 1.e-15);
    EXPECT_NEAR(quadrature[5][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[6][0], -.4, 1.e-15);
    EXPECT_NEAR(quadrature[6][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[7][0], -.3, 1.e-15);
    EXPECT_NEAR(quadrature[7][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[8][0], -.2, 1.e-15);
    EXPECT_NEAR(quadrature[8][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[9][0], -.1, 1.e-15);
    EXPECT_NEAR(quadrature[9][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[10][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[10][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[11][0], .1, 1.e-15);
    EXPECT_NEAR(quadrature[11][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[12][0], .2, 1.e-15);
    EXPECT_NEAR(quadrature[12][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[13][0], .3, 1.e-15);
    EXPECT_NEAR(quadrature[13][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[14][0], .4, 1.e-15);
    EXPECT_NEAR(quadrature[14][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[15][0], .5, 1.e-15);
    EXPECT_NEAR(quadrature[15][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[16][0], .6, 1.e-15);
    EXPECT_NEAR(quadrature[16][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[17][0], .7, 1.e-15);
    EXPECT_NEAR(quadrature[17][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[18][0], .8, 1.e-15);
    EXPECT_NEAR(quadrature[18][1], (2. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[19][0], .9, 1.e-15);
    EXPECT_NEAR(quadrature[19][1], (4. / 3.) / 10., 1.e-15);

    EXPECT_NEAR(quadrature[20][0], 1., 1.e-15);
    EXPECT_NEAR(quadrature[20][1], (1. / 3.) / 10., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GL_TwoSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.quadrature_rule =
        kynema::interfaces::components::BeamInput::QuadratureRule::GaussLegendre;
    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], 0., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 2., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GL_ThreeSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.quadrature_rule =
        kynema::interfaces::components::BeamInput::QuadratureRule::GaussLegendre;
    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -.5, 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], .5, 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 1., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GL_FourSections_NoRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.75, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.quadrature_rule =
        kynema::interfaces::components::BeamInput::QuadratureRule::GaussLegendre;
    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -.5, 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], .25, 1.e-15);
    EXPECT_NEAR(quadrature[1][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[2][0], .75, 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .5, 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GL_TwoSections_OneRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.quadrature_rule =
        kynema::interfaces::components::BeamInput::QuadratureRule::GaussLegendre;
    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -1. / std::sqrt(3.), 1.e-15);
    EXPECT_NEAR(quadrature[0][1], 1., 1.e-15);

    EXPECT_NEAR(quadrature[1][0], 1. / std::sqrt(3.), 1.e-15);
    EXPECT_NEAR(quadrature[1][1], 1., 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GL_ThreeSections_OneRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.quadrature_rule =
        kynema::interfaces::components::BeamInput::QuadratureRule::GaussLegendre;
    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -.5 - 1. / std::sqrt(3.) / 2., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -.5 + 1. / std::sqrt(3.) / 2., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[2][0], .5 - 1. / std::sqrt(3.) / 2., 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .5, 1.e-15);

    EXPECT_NEAR(quadrature[3][0], .5 + 1. / std::sqrt(3.) / 2., 1.e-15);
    EXPECT_NEAR(quadrature[3][1], .5, 1.e-15);
}

TEST(BeamComponentTest, Quadrature_GL_ElevenSections_OneRefinement) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    auto mass_stiff_array =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.1, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.2, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.3, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.4, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.5, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.6, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.7, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.8, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(.9, mass_stiff_array, mass_stiff_array),
        kynema::interfaces::components::Section(1., mass_stiff_array, mass_stiff_array)
    };

    beam_input.quadrature_rule =
        kynema::interfaces::components::BeamInput::QuadratureRule::GaussLegendre;
    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto quadrature = beam_elements.front().quadrature;

    EXPECT_NEAR(quadrature[0][0], -.9 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[0][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[1][0], -.9 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[1][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[2][0], -.7 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[2][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[3][0], -.7 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[3][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[4][0], -.5 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[4][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[5][0], -.5 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[5][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[6][0], -.3 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[6][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[7][0], -.3 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[7][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[8][0], -.1 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[8][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[9][0], -.1 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[9][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[10][0], .1 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[10][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[11][0], .1 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[11][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[12][0], .3 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[12][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[13][0], .3 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[13][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[14][0], .5 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[14][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[15][0], .5 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[15][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[16][0], .7 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[16][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[17][0], .7 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[17][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[18][0], .9 - 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[18][1], .1, 1.e-15);

    EXPECT_NEAR(quadrature[19][0], .9 + 1. / std::sqrt(3.) / 10., 1.e-15);
    EXPECT_NEAR(quadrature[19][1], .1, 1.e-15);
}
}  // namespace kynema::tests
