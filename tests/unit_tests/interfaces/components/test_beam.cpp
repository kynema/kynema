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

TEST(BeamComponentTest, Section_GLL_TwoSections_NoRefinement_NoTwist_ConstantSection) {
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
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 2UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, 1., 1.e-15);

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array[i][j], 1.e-15);
        }
    }
}

TEST(BeamComponentTest, Section_GLL_TwoSections_NoRefinement_NoTwist_LinearSection) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    constexpr auto mass_stiff_array_0 =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    constexpr auto mass_stiff_array_1 =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 2., 0., 0., 0.}, std::array{0., 0., 0., 2., 0., 0.},
                   std::array{0., 0., 0., 0., 2., 0.}, std::array{0., 0., 0., 0., 0., 2.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array_0, mass_stiff_array_0),
        kynema::interfaces::components::Section(1., mass_stiff_array_1, mass_stiff_array_1)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 2UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, 1., 1.e-15);

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array_0[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array_1[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array_0[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array_1[i][j], 1.e-15);
        }
    }
}

TEST(BeamComponentTest, Section_GLL_ThreeSections_NoRefinement_NoTwist_LinearSection) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    constexpr auto mass_stiff_array_0 =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    constexpr auto mass_stiff_array_1 =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 2., 0., 0., 0.}, std::array{0., 0., 0., 2., 0., 0.},
                   std::array{0., 0., 0., 0., 2., 0.}, std::array{0., 0., 0., 0., 0., 2.}};
    constexpr auto mass_stiff_array_2 =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 4., 0., 0., 0., 0.},
                   std::array{0., 0., 4., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 4., 0.}, std::array{0., 0., 0., 0., 0., 4.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array_0, mass_stiff_array_0),
        kynema::interfaces::components::Section(.5, mass_stiff_array_1, mass_stiff_array_1),
        kynema::interfaces::components::Section(1., mass_stiff_array_2, mass_stiff_array_2)
    };

    beam_input.section_refinement = 0UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 3UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, .5, 1.e-15);
    EXPECT_NEAR(sections[2].position, 1., 1.e-15);

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array_0[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array_1[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].M_star[i][j], mass_stiff_array_2[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array_0[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array_1[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].C_star[i][j], mass_stiff_array_2[i][j], 1.e-15);
        }
    }
}

TEST(BeamComponentTest, Section_GLL_ThreeSections_OneRefinement_NoTwist_LinearSection) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 0.};

    constexpr auto mass_stiff_array_0 =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., 1.}};
    constexpr auto mass_stiff_array_1 =
        std::array{std::array{1.5, 0., 0., 0., 0., 0.}, std::array{0., 1.5, 0., 0., 0., 0.},
                   std::array{0., 0., 1.5, 0., 0., 0.}, std::array{0., 0., 0., 1.5, 0., 0.},
                   std::array{0., 0., 0., 0., 1.5, 0.}, std::array{0., 0., 0., 0., 0., 1.5}};
    constexpr auto mass_stiff_array_2 =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 2., 0., 0., 0.}, std::array{0., 0., 0., 2., 0., 0.},
                   std::array{0., 0., 0., 0., 2., 0.}, std::array{0., 0., 0., 0., 0., 2.}};
    constexpr auto mass_stiff_array_3 =
        std::array{std::array{3., 0., 0., 0., 0., 0.}, std::array{0., 3., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 3., 0., 0.},
                   std::array{0., 0., 0., 0., 3., 0.}, std::array{0., 0., 0., 0., 0., 3.}};
    constexpr auto mass_stiff_array_4 =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 4., 0., 0., 0., 0.},
                   std::array{0., 0., 4., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 4., 0.}, std::array{0., 0., 0., 0., 0., 4.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array_0, mass_stiff_array_0),
        kynema::interfaces::components::Section(.5, mass_stiff_array_2, mass_stiff_array_2),
        kynema::interfaces::components::Section(1., mass_stiff_array_4, mass_stiff_array_4)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 5UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, .25, 1.e-15);
    EXPECT_NEAR(sections[2].position, .5, 1.e-15);
    EXPECT_NEAR(sections[3].position, .75, 1.e-15);
    EXPECT_NEAR(sections[4].position, 1., 1.e-15);

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array_0[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array_1[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].M_star[i][j], mass_stiff_array_2[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].M_star[i][j], mass_stiff_array_3[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].M_star[i][j], mass_stiff_array_4[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array_0[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array_1[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].C_star[i][j], mass_stiff_array_2[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].C_star[i][j], mass_stiff_array_3[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].C_star[i][j], mass_stiff_array_4[i][j], 1.e-15);
        }
    }
}

TEST(BeamComponentTest, Section_GLL_ThreeSections_OneRefinement_ConstantTwist_ConstantSection) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{45., 45.};

    constexpr auto mass_stiff_array_0 =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 5., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_2 =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 3., 0., 0., 0., 0.},
                   std::array{0., 0., 4., 0., 0., 0.}, std::array{0., 0., 0., 5., 0., 0.},
                   std::array{0., 0., 0., 0., 6., 0.}, std::array{0., 0., 0., 0., 0., 7.}};
    constexpr auto mass_stiff_array_4 =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 5., 0., 0., 0., 0.},
                   std::array{0., 0., 6., 0., 0., 0.}, std::array{0., 0., 0., 7., 0., 0.},
                   std::array{0., 0., 0., 0., 8., 0.}, std::array{0., 0., 0., 0., 0., 9.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array_0, mass_stiff_array_0),
        kynema::interfaces::components::Section(.5, mass_stiff_array_2, mass_stiff_array_2),
        kynema::interfaces::components::Section(1., mass_stiff_array_4, mass_stiff_array_4)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 5UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, .25, 1.e-15);
    EXPECT_NEAR(sections[2].position, .5, 1.e-15);
    EXPECT_NEAR(sections[3].position, .75, 1.e-15);
    EXPECT_NEAR(sections[4].position, 1., 1.e-15);

    constexpr auto mass_stiff_array_0_exact =
        std::array{std::array{1., 0., 0., 0., 0., 0.},   std::array{0., 2.5, -.5, 0., 0., 0.},
                   std::array{0., -.5, 2.5, 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 5.5, -.5}, std::array{0., 0., 0., 0., -.5, 5.5}};
    constexpr auto mass_stiff_array_1_exact =
        std::array{std::array{1.5, 0., 0., 0., 0., 0.}, std::array{0., 3., -.5, 0., 0., 0.},
                   std::array{0., -.5, 3., 0., 0., 0.}, std::array{0., 0., 0., 4.5, 0., 0.},
                   std::array{0., 0., 0., 0., 6., -.5}, std::array{0., 0., 0., 0., -.5, 6.}};
    constexpr auto mass_stiff_array_2_exact =
        std::array{std::array{2., 0., 0., 0., 0., 0.},   std::array{0., 3.5, -.5, 0., 0., 0.},
                   std::array{0., -.5, 3.5, 0., 0., 0.}, std::array{0., 0., 0., 5., 0., 0.},
                   std::array{0., 0., 0., 0., 6.5, -.5}, std::array{0., 0., 0., 0., -.5, 6.5}};
    constexpr auto mass_stiff_array_3_exact =
        std::array{std::array{3., 0., 0., 0., 0., 0.},   std::array{0., 4.5, -.5, 0., 0., 0.},
                   std::array{0., -.5, 4.5, 0., 0., 0.}, std::array{0., 0., 0., 6., 0., 0.},
                   std::array{0., 0., 0., 0., 7.5, -.5}, std::array{0., 0., 0., 0., -.5, 7.5}};
    constexpr auto mass_stiff_array_4_exact =
        std::array{std::array{4., 0., 0., 0., 0., 0.},   std::array{0., 5.5, -.5, 0., 0., 0.},
                   std::array{0., -.5, 5.5, 0., 0., 0.}, std::array{0., 0., 0., 7., 0., 0.},
                   std::array{0., 0., 0., 0., 8.5, -.5}, std::array{0., 0., 0., 0., -.5, 8.5}};

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array_0_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array_1_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].M_star[i][j], mass_stiff_array_2_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].M_star[i][j], mass_stiff_array_3_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].M_star[i][j], mass_stiff_array_4_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array_0_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array_1_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].C_star[i][j], mass_stiff_array_2_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].C_star[i][j], mass_stiff_array_3_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].C_star[i][j], mass_stiff_array_4_exact[i][j], 1.e-15);
        }
    }
}

TEST(BeamComponentTest, Section_GLL_ThreeSections_OneRefinement_LinearTwist_ConstantSection) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., 1.0};
    beam_input.ref_axis.twist = std::vector{0., 180.};

    constexpr auto mass_stiff_array_0 =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 5., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_2 =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 3., 0., 0., 0., 0.},
                   std::array{0., 0., 4., 0., 0., 0.}, std::array{0., 0., 0., 5., 0., 0.},
                   std::array{0., 0., 0., 0., 6., 0.}, std::array{0., 0., 0., 0., 0., 7.}};
    constexpr auto mass_stiff_array_4 =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 5., 0., 0., 0., 0.},
                   std::array{0., 0., 6., 0., 0., 0.}, std::array{0., 0., 0., 7., 0., 0.},
                   std::array{0., 0., 0., 0., 8., 0.}, std::array{0., 0., 0., 0., 0., 9.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array_0, mass_stiff_array_0),
        kynema::interfaces::components::Section(.5, mass_stiff_array_2, mass_stiff_array_2),
        kynema::interfaces::components::Section(1., mass_stiff_array_4, mass_stiff_array_4)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 5UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, .25, 1.e-15);
    EXPECT_NEAR(sections[2].position, .5, 1.e-15);
    EXPECT_NEAR(sections[3].position, .75, 1.e-15);
    EXPECT_NEAR(sections[4].position, 1., 1.e-15);

    constexpr auto mass_stiff_array_0_exact =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 5., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_1_exact =
        std::array{std::array{1.5, 0., 0., 0., 0., 0.}, std::array{0., 3., -.5, 0., 0., 0.},
                   std::array{0., -.5, 3., 0., 0., 0.}, std::array{0., 0., 0., 4.5, 0., 0.},
                   std::array{0., 0., 0., 0., 6., -.5}, std::array{0., 0., 0., 0., -.5, 6.}};
    constexpr auto mass_stiff_array_2_exact =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 4., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 5., 0., 0.},
                   std::array{0., 0., 0., 0., 7., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_3_exact =
        std::array{std::array{3., 0., 0., 0., 0., 0.},  std::array{0., 4.5, .5, 0., 0., 0.},
                   std::array{0., .5, 4.5, 0., 0., 0.}, std::array{0., 0., 0., 6., 0., 0.},
                   std::array{0., 0., 0., 0., 7.5, .5}, std::array{0., 0., 0., 0., .5, 7.5}};
    constexpr auto mass_stiff_array_4_exact =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 5., 0., 0., 0., 0.},
                   std::array{0., 0., 6., 0., 0., 0.}, std::array{0., 0., 0., 7., 0., 0.},
                   std::array{0., 0., 0., 0., 8., 0.}, std::array{0., 0., 0., 0., 0., 9.}};

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array_0_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array_1_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].M_star[i][j], mass_stiff_array_2_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].M_star[i][j], mass_stiff_array_3_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].M_star[i][j], mass_stiff_array_4_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array_0_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array_1_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].C_star[i][j], mass_stiff_array_2_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].C_star[i][j], mass_stiff_array_3_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].C_star[i][j], mass_stiff_array_4_exact[i][j], 1.e-15);
        }
    }
}

TEST(BeamComponentTest, Section_GLL_ThreeSections_OneRefinement_BiLinearTwist_ConstantSection) {
    auto model = kynema::Model();
    auto beam_input = kynema::interfaces::components::BeamInput{};

    beam_input.element_order = 2UL;

    beam_input.ref_axis.coordinate_grid = std::vector{0., 1.0};
    beam_input.ref_axis.coordinates = std::vector{std::array{0., 0., 0.}, std::array{1., 0., 0.}};
    beam_input.ref_axis.twist_grid = std::vector{0., .5, 1.0};
    beam_input.ref_axis.twist = std::vector{0., 90., 0.};

    constexpr auto mass_stiff_array_0 =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 5., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_2 =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 3., 0., 0., 0., 0.},
                   std::array{0., 0., 4., 0., 0., 0.}, std::array{0., 0., 0., 5., 0., 0.},
                   std::array{0., 0., 0., 0., 6., 0.}, std::array{0., 0., 0., 0., 0., 7.}};
    constexpr auto mass_stiff_array_4 =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 5., 0., 0., 0., 0.},
                   std::array{0., 0., 6., 0., 0., 0.}, std::array{0., 0., 0., 7., 0., 0.},
                   std::array{0., 0., 0., 0., 8., 0.}, std::array{0., 0., 0., 0., 0., 9.}};
    beam_input.sections = std::vector{
        kynema::interfaces::components::Section(0., mass_stiff_array_0, mass_stiff_array_0),
        kynema::interfaces::components::Section(.5, mass_stiff_array_2, mass_stiff_array_2),
        kynema::interfaces::components::Section(1., mass_stiff_array_4, mass_stiff_array_4)
    };

    beam_input.section_refinement = 1UL;

    const auto beam = kynema::interfaces::components::Beam(beam_input, model);
    const auto beam_elements = model.GetBeamElements();
    const auto& sections = beam_elements.front().sections;

    EXPECT_EQ(sections.size(), 5UL);

    EXPECT_NEAR(sections[0].position, 0., 1.e-15);
    EXPECT_NEAR(sections[1].position, .25, 1.e-15);
    EXPECT_NEAR(sections[2].position, .5, 1.e-15);
    EXPECT_NEAR(sections[3].position, .75, 1.e-15);
    EXPECT_NEAR(sections[4].position, 1., 1.e-15);

    constexpr auto mass_stiff_array_0_exact =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 2., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 4., 0., 0.},
                   std::array{0., 0., 0., 0., 5., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_1_exact =
        std::array{std::array{1.5, 0., 0., 0., 0., 0.}, std::array{0., 3., -.5, 0., 0., 0.},
                   std::array{0., -.5, 3., 0., 0., 0.}, std::array{0., 0., 0., 4.5, 0., 0.},
                   std::array{0., 0., 0., 0., 6., -.5}, std::array{0., 0., 0., 0., -.5, 6.}};
    constexpr auto mass_stiff_array_2_exact =
        std::array{std::array{2., 0., 0., 0., 0., 0.}, std::array{0., 4., 0., 0., 0., 0.},
                   std::array{0., 0., 3., 0., 0., 0.}, std::array{0., 0., 0., 5., 0., 0.},
                   std::array{0., 0., 0., 0., 7., 0.}, std::array{0., 0., 0., 0., 0., 6.}};
    constexpr auto mass_stiff_array_3_exact =
        std::array{std::array{3., 0., 0., 0., 0., 0.},   std::array{0., 4.5, -.5, 0., 0., 0.},
                   std::array{0., -.5, 4.5, 0., 0., 0.}, std::array{0., 0., 0., 6., 0., 0.},
                   std::array{0., 0., 0., 0., 7.5, -.5}, std::array{0., 0., 0., 0., -.5, 7.5}};
    constexpr auto mass_stiff_array_4_exact =
        std::array{std::array{4., 0., 0., 0., 0., 0.}, std::array{0., 5., 0., 0., 0., 0.},
                   std::array{0., 0., 6., 0., 0., 0.}, std::array{0., 0., 0., 7., 0., 0.},
                   std::array{0., 0., 0., 0., 8., 0.}, std::array{0., 0., 0., 0., 0., 9.}};

    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(sections[0].M_star[i][j], mass_stiff_array_0_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].M_star[i][j], mass_stiff_array_1_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].M_star[i][j], mass_stiff_array_2_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].M_star[i][j], mass_stiff_array_3_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].M_star[i][j], mass_stiff_array_4_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[0].C_star[i][j], mass_stiff_array_0_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[1].C_star[i][j], mass_stiff_array_1_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[2].C_star[i][j], mass_stiff_array_2_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[3].C_star[i][j], mass_stiff_array_3_exact[i][j], 1.e-15);
            EXPECT_NEAR(sections[4].C_star[i][j], mass_stiff_array_4_exact[i][j], 1.e-15);
        }
    }
}
}  // namespace kynema::tests
