#include <vector>

#include <gtest/gtest.h>

#include "interfaces/components/aerodynamics.hpp"
#include "math/interpolation.hpp"

namespace kynema::tests {

TEST(AerodynamicsComponentTest, CalculateAngleOfAttack_Case1) {
    constexpr auto flow_angle = 0.;
    constexpr auto expected_aoa = 0.;
    const auto v_rel = std::array{0., std::cos(flow_angle), std::sin(flow_angle)};
    const auto aoa = kynema::interfaces::components::CalculateAngleOfAttack(v_rel);

    EXPECT_NEAR(aoa, expected_aoa, 1.e-15);
}

TEST(AerodynamicsComponentTest, CalculateAngleOfAttack_Case2) {
    constexpr auto flow_angle = -.1;
    constexpr auto expected_aoa = .1;
    const auto v_rel = std::array{0., std::cos(flow_angle), std::sin(flow_angle)};
    const auto aoa = kynema::interfaces::components::CalculateAngleOfAttack(v_rel);

    EXPECT_NEAR(aoa, expected_aoa, 1.e-15);
}

TEST(AerodynamicsComponentTest, CalculateAngleOfAttack_Case3) {
    constexpr auto flow_angle = .2;
    constexpr auto expected_aoa = -.2;
    const auto v_rel = std::array{0., std::cos(flow_angle), std::sin(flow_angle)};
    const auto aoa = kynema::interfaces::components::CalculateAngleOfAttack(v_rel);

    EXPECT_NEAR(aoa, expected_aoa, 1.e-15);
}

TEST(AerodynamicsComponentTest, CalculateAngleOfAttack_Case4) {
    constexpr auto flow_angle = 1. - M_PI;
    constexpr auto expected_aoa = M_PI - 1.;
    const auto v_rel = std::array{0., std::cos(flow_angle), std::sin(flow_angle)};
    const auto aoa = kynema::interfaces::components::CalculateAngleOfAttack(v_rel);

    EXPECT_NEAR(aoa, expected_aoa, 1.e-15);
}

TEST(AerodynamicsComponentTest, CalculateAerodynamicLoad_Case1) {
    const auto aoa_polar = std::vector{-1., 1.};
    const auto cl_polar = std::vector{0., 1.};
    const auto cd_polar = std::vector{.5, 0.};
    const auto cm_polar = std::vector{-0.01, -0.03};

    constexpr auto v_inflow = std::array{0., 10., 0.};
    constexpr auto v_motion = std::array{0., 0., 0.};
    constexpr auto chord = 2.;
    constexpr auto delta_s = 1.5;
    constexpr auto fluid_density = 1.225;
    constexpr auto qqr = std::array{1., 0., 0., 0.};
    constexpr auto con_force = std::array{0., -.5, 0.};

    double alpha, cn, ct, cm;
    std::array<double, 3> v_rel;

    auto ref_axis_moment = std::array<double, 3>{};
    auto load = kynema::interfaces::components::CalculateAerodynamicLoad(
        ref_axis_moment, v_inflow, v_motion, aoa_polar, cl_polar, cd_polar, cm_polar, chord, delta_s,
        fluid_density, con_force, qqr, v_rel, alpha, cn, ct, cm
    );

    constexpr auto expected_load = std::array{0., 45.9375, -91.875, -7.35, 0., 0.};
    constexpr auto expected_ref_axis_moment = std::array{-7.35 - 91.875 / 2., 0., 0.};

    EXPECT_NEAR(load[0], expected_load[0], 1.e-10);
    EXPECT_NEAR(load[1], expected_load[1], 1.e-10);
    EXPECT_NEAR(load[2], expected_load[2], 1.e-10);
    EXPECT_NEAR(load[3], expected_load[3], 1.e-10);
    EXPECT_NEAR(load[4], expected_load[4], 1.e-10);
    EXPECT_NEAR(load[5], expected_load[5], 1.e-10);

    EXPECT_NEAR(ref_axis_moment[0], expected_ref_axis_moment[0], 1.e-10);
    EXPECT_NEAR(ref_axis_moment[1], expected_ref_axis_moment[1], 1.e-10);
    EXPECT_NEAR(ref_axis_moment[2], expected_ref_axis_moment[2], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateAerodynamicLoad_Case2) {
    const auto aoa_polar = std::vector{-1., 1.};
    const auto cl_polar = std::vector{0., 1.};
    const auto cd_polar = std::vector{.5, 0.};
    const auto cm_polar = std::vector{-0.01, -0.03};

    constexpr auto v_inflow = std::array{0., 9.950041652780259, -0.9983341664682815};
    constexpr auto v_motion = std::array{0., 0., 0.};
    constexpr auto chord = 2.;
    constexpr auto delta_s = 1.5;
    constexpr auto fluid_density = 1.225;
    constexpr auto qqr = std::array{1., 0., 0., 0.};
    constexpr auto con_force = std::array{0., -.5, 0.};

    double alpha, cn, ct, cm;
    std::array<double, 3> v_rel;

    auto ref_axis_moment = std::array<double, 3>{};
    auto load = kynema::interfaces::components::CalculateAerodynamicLoad(
        ref_axis_moment, v_inflow, v_motion, aoa_polar, cl_polar, cd_polar, cm_polar, chord, delta_s,
        fluid_density, con_force, qqr, v_rel, alpha, cn, ct, cm
    );

    constexpr auto expected_load =
        std::array{0., 31.04778878834331, -104.68509627290281, -7.7175, 0., 0.};
    constexpr auto expected_ref_axis_moment = std::array{-7.7175 - 104.68509627290281 / 2., 0., 0.};

    EXPECT_NEAR(load[0], expected_load[0], 1.e-10);
    EXPECT_NEAR(load[1], expected_load[1], 1.e-10);
    EXPECT_NEAR(load[2], expected_load[2], 1.e-10);
    EXPECT_NEAR(load[3], expected_load[3], 1.e-10);
    EXPECT_NEAR(load[4], expected_load[4], 1.e-10);
    EXPECT_NEAR(load[5], expected_load[5], 1.e-10);

    EXPECT_NEAR(ref_axis_moment[0], expected_ref_axis_moment[0], 1.e-10);
    EXPECT_NEAR(ref_axis_moment[1], expected_ref_axis_moment[1], 1.e-10);
    EXPECT_NEAR(ref_axis_moment[2], expected_ref_axis_moment[2], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateConMotionVector_Case1) {
    constexpr auto ac_to_ref_axis_horizontal = 1.;
    constexpr auto chord_to_ref_axis_vertical = 0.;

    const auto ac_vector = kynema::interfaces::components::CalculateConMotionVector(
        ac_to_ref_axis_horizontal, chord_to_ref_axis_vertical
    );

    constexpr auto expected_ac_vector = std::array{0., -1., 0.};
    EXPECT_NEAR(ac_vector[0], expected_ac_vector[0], 1.e-10);
    EXPECT_NEAR(ac_vector[1], expected_ac_vector[1], 1.e-10);
    EXPECT_NEAR(ac_vector[2], expected_ac_vector[2], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateConMotionVector_Case2) {
    constexpr auto ac_to_ref_axis_horizontal = 1.;
    constexpr auto chord_to_ref_axis_vertical = .5;

    const auto ac_vector = kynema::interfaces::components::CalculateConMotionVector(
        ac_to_ref_axis_horizontal, chord_to_ref_axis_vertical
    );

    constexpr auto expected_ac_vector = std::array{0., -1., .5};
    EXPECT_NEAR(ac_vector[0], expected_ac_vector[0], 1.e-10);
    EXPECT_NEAR(ac_vector[1], expected_ac_vector[1], 1.e-10);
    EXPECT_NEAR(ac_vector[2], expected_ac_vector[2], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateJoacibanXi_Case1) {
    const auto aero_node_xi = std::vector{-1., 0., 1.};
    const auto jacobian_xi = kynema::interfaces::components::CalculateJacobianXi(aero_node_xi);

    constexpr auto expected_jacobian_xi = std::array{-1., -.75, -.5, 0., .5, .75, 1.};
    EXPECT_NEAR(jacobian_xi[0], expected_jacobian_xi[0], 1.e-10);
    EXPECT_NEAR(jacobian_xi[1], expected_jacobian_xi[1], 1.e-10);
    EXPECT_NEAR(jacobian_xi[2], expected_jacobian_xi[2], 1.e-10);
    EXPECT_NEAR(jacobian_xi[3], expected_jacobian_xi[3], 1.e-10);
    EXPECT_NEAR(jacobian_xi[4], expected_jacobian_xi[4], 1.e-10);
    EXPECT_NEAR(jacobian_xi[5], expected_jacobian_xi[5], 1.e-10);
    EXPECT_NEAR(jacobian_xi[6], expected_jacobian_xi[6], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateJoacibanXi_Case2) {
    const auto aero_node_xi = std::vector{-1., .2, 1.};
    const auto jacobian_xi = kynema::interfaces::components::CalculateJacobianXi(aero_node_xi);

    constexpr auto expected_jacobian_xi = std::array{-1., -.7, -.4, .2, .6, .8, 1.};
    EXPECT_NEAR(jacobian_xi[0], expected_jacobian_xi[0], 1.e-10);
    EXPECT_NEAR(jacobian_xi[1], expected_jacobian_xi[1], 1.e-10);
    EXPECT_NEAR(jacobian_xi[2], expected_jacobian_xi[2], 1.e-10);
    EXPECT_NEAR(jacobian_xi[3], expected_jacobian_xi[3], 1.e-10);
    EXPECT_NEAR(jacobian_xi[4], expected_jacobian_xi[4], 1.e-10);
    EXPECT_NEAR(jacobian_xi[5], expected_jacobian_xi[5], 1.e-10);
    EXPECT_NEAR(jacobian_xi[6], expected_jacobian_xi[6], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateJoacibanXi_Case3) {
    const auto aero_node_xi = std::vector{-1., -.2, .2, 1.};
    const auto jacobian_xi = kynema::interfaces::components::CalculateJacobianXi(aero_node_xi);

    constexpr auto expected_jacobian_xi = std::array{-1., -.8, -.6, -.2, 0., .2, .6, .8, 1.};
    EXPECT_NEAR(jacobian_xi[0], expected_jacobian_xi[0], 1.e-10);
    EXPECT_NEAR(jacobian_xi[1], expected_jacobian_xi[1], 1.e-10);
    EXPECT_NEAR(jacobian_xi[2], expected_jacobian_xi[2], 1.e-10);
    EXPECT_NEAR(jacobian_xi[3], expected_jacobian_xi[3], 1.e-10);
    EXPECT_NEAR(jacobian_xi[4], expected_jacobian_xi[4], 1.e-10);
    EXPECT_NEAR(jacobian_xi[5], expected_jacobian_xi[5], 1.e-10);
    EXPECT_NEAR(jacobian_xi[6], expected_jacobian_xi[6], 1.e-10);
    EXPECT_NEAR(jacobian_xi[7], expected_jacobian_xi[7], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateAeroNodeWidths_Straight) {
    const auto beam_node_xi = std::vector{-1., -.75, -.5, -.25, .0, .25, .5, .75, 1.};
    const auto node_x = std::vector{0., 1., 2., 3., 4., 5., 6., 7., 8., 0., 0., 0., 0., 0.,
                                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    const auto aero_node_xi = std::vector{-1., 0., 1.};
    const auto jacobian_xi = kynema::interfaces::components::CalculateJacobianXi(aero_node_xi);

    auto jacobian_integration_matrix = std::vector<double>(beam_node_xi.size() * jacobian_xi.size());

    auto weights = std::vector<double>{};
    for (auto i = 0U; i < jacobian_xi.size(); ++i) {
        kynema::math::LagrangePolynomialDerivWeights(jacobian_xi[i], beam_node_xi, weights);
        for (auto j = 0U; j < beam_node_xi.size(); ++j) {
            jacobian_integration_matrix[i * beam_node_xi.size() + j] = weights[j];
        }
    }

    const auto widths = kynema::interfaces::components::CalculateAeroNodeWidths(
        jacobian_xi, jacobian_integration_matrix, node_x
    );

    EXPECT_EQ(widths.size(), 3UL);

    constexpr auto expected_widths = std::array{.25 * 8., .5 * 8., .25 * 8.};
    EXPECT_NEAR(widths[0], expected_widths[0], 1.e-10);
    EXPECT_NEAR(widths[1], expected_widths[1], 1.e-10);
    EXPECT_NEAR(widths[2], expected_widths[2], 1.e-10);
}

TEST(AerodynamicsComponentTest, CalculateAeroNodeWidths_Curved) {
    const auto beam_node_xi = std::vector{-1., 0., 1.};
    const auto node_x = std::vector{0., 1., 2., 0., 1., .5, 0., 0., 0.};

    const auto aero_node_xi = std::vector{-1., -.5, 0., .5, 1.};
    const auto jacobian_xi = kynema::interfaces::components::CalculateJacobianXi(aero_node_xi);

    auto jacobian_integration_matrix = std::vector<double>(beam_node_xi.size() * jacobian_xi.size());

    auto weights = std::vector<double>{};
    for (auto i = 0U; i < jacobian_xi.size(); ++i) {
        kynema::math::LagrangePolynomialDerivWeights(jacobian_xi[i], beam_node_xi, weights);
        for (auto j = 0U; j < beam_node_xi.size(); ++j) {
            jacobian_integration_matrix[i * beam_node_xi.size() + j] = weights[j];
        }
    }

    const auto widths = kynema::interfaces::components::CalculateAeroNodeWidths(
        jacobian_xi, jacobian_integration_matrix, node_x
    );

    EXPECT_EQ(widths.size(), 5UL);

    constexpr auto expected_widths = std::array{
        .46400663710393675, .71135714204928224, .52584462380517116, .56739053334406075,
        .36524408545476744
    };
    EXPECT_NEAR(widths[0], expected_widths[0], 1.e-10);
    EXPECT_NEAR(widths[1], expected_widths[1], 1.e-10);
    EXPECT_NEAR(widths[2], expected_widths[2], 1.e-10);
    EXPECT_NEAR(widths[3], expected_widths[3], 1.e-10);
    EXPECT_NEAR(widths[4], expected_widths[4], 1.e-10);
}
}  // namespace kynema::tests
