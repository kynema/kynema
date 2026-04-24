#include <array>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>

#include <gtest/gtest.h>

#include "elements/beams/beam_section.hpp"
#include "elements/beams/hollow_circle_properties.hpp"

namespace kynema_fmb::beams::tests {

class HollowCirclePropertiesTest : public ::testing::Test {
protected:
    static constexpr double kTolerance = 1e-12;

    struct TestParameters {
        // Geometric properties
        double outer_diameter{1.};
        double wall_thickness{0.1};
        double nu{0.33};

        // Material properties
        double E{200.E9};   // Steel Young's modulus (Pa)
        double G{80.E9};    // Steel shear modulus (Pa)
        double rho{7850.};  // Steel density (kg/m³)

        // Offset parameters
        double x_C{0.};
        double y_C{0.};
        double theta_p{0.};
        double x_S{0.};
        double y_S{0.};
        double theta_s{0.};
        double x_G{0.};
        double y_G{0.};
        double theta_i{0.};

        // Position parameter
        double s{0.5};
    };

    TestParameters tp_;

    static void ExpectMatrixEqual(
        const std::array<std::array<double, 6>, 6>& actual,
        const std::array<std::array<double, 6>, 6>& expected
    ) {
        for (auto i : std::views::iota(0U, 6U)) {
            for (auto j : std::views::iota(0U, 6U)) {
                EXPECT_NEAR(actual[i][j], expected[i][j], kTolerance)
                    << "Matrix element [" << i << "][" << j << "] differs";
            }
        }
    }
};

TEST_F(HollowCirclePropertiesTest, InvalidGeometryThrows) {
    // Test that invalid geometry throws an exception
    EXPECT_THROW(
        CalculateHollowCircleProperties(1., 0.6, 0.3),  // thickness = radius
        std::invalid_argument
    );

    EXPECT_THROW(
        CalculateHollowCircleProperties(1., 0.5, 0.3),  // thickness > radius
        std::invalid_argument
    );
}

TEST_F(HollowCirclePropertiesTest, CalculateGeometricProperties) {
    // Test basic geometric property calculations
    auto properties =
        CalculateHollowCircleProperties(tp_.outer_diameter, tp_.wall_thickness, tp_.nu);

    // Calculate expected values manually
    const double outer_radius = tp_.outer_diameter / 2.;
    const double inner_radius = outer_radius - tp_.wall_thickness;

    const double expected_area =
        std::numbers::pi * (std::pow(outer_radius, 2) - std::pow(inner_radius, 2));
    const double expected_Ixx =
        std::numbers::pi * (std::pow(outer_radius, 4) - std::pow(inner_radius, 4)) / 4.;
    const double expected_J =
        std::numbers::pi * (std::pow(outer_radius, 4) - std::pow(inner_radius, 4)) / 2.;
    const double expected_kx = (6. * (1. + tp_.nu)) / (7. + 6. * tp_.nu);

    EXPECT_NEAR(properties.area, expected_area, kTolerance);
    EXPECT_NEAR(properties.Ixx, expected_Ixx, kTolerance);
    EXPECT_NEAR(properties.Iyy, expected_Ixx, kTolerance);  // Circular symmetry
    EXPECT_NEAR(properties.J, expected_J, kTolerance);
    EXPECT_NEAR(properties.kx, expected_kx, kTolerance);
    EXPECT_NEAR(properties.ky, expected_kx, kTolerance);  // Circular symmetry
}

TEST_F(HollowCirclePropertiesTest, SectionalMatrices_Uncoupled) {
    // Test stiffness and mass matrix generation with no coupling
    auto beam_section = GenerateHollowCircleSection(
        tp_.s, tp_.E, tp_.G, tp_.rho, tp_.outer_diameter, tp_.wall_thickness, tp_.nu, tp_.x_C,
        tp_.y_C, tp_.theta_p, tp_.x_S, tp_.y_S, tp_.theta_s, tp_.x_G, tp_.y_G, tp_.theta_i
    );

    auto properties =
        CalculateHollowCircleProperties(tp_.outer_diameter, tp_.wall_thickness, tp_.nu);

    // Calculate expected stiffness values
    const double EA = tp_.E * properties.area;
    const double EI_x = tp_.E * properties.Ixx;
    const double EI_y = tp_.E * properties.Iyy;
    const double GKt = tp_.G * properties.J;
    const double GA = tp_.G * properties.area;

    auto expected_stiffness = std::array<std::array<double, 6>, 6>{};
    expected_stiffness[0][0] = GA * properties.kx;  // Shear stiffness in x
    expected_stiffness[1][1] = GA * properties.ky;  // Shear stiffness in y
    expected_stiffness[2][2] = EA;                  // Axial stiffness
    expected_stiffness[3][3] = EI_x;                // Bending stiffness about x
    expected_stiffness[4][4] = EI_y;                // Bending stiffness about y
    expected_stiffness[5][5] = GKt;                 // Torsional stiffness

    ExpectMatrixEqual(beam_section.C_star, expected_stiffness);

    // Calculate expected mass values
    const double m = tp_.rho * properties.area;
    const double I_x = tp_.rho * properties.Ixx;
    const double I_y = tp_.rho * properties.Iyy;
    const double I_p = I_x + I_y;

    auto expected_mass = std::array<std::array<double, 6>, 6>{};
    expected_mass[0][0] = m;    // Translational mass in x
    expected_mass[1][1] = m;    // Translational mass in y
    expected_mass[2][2] = m;    // Translational mass in z
    expected_mass[3][3] = I_x;  // Rotational inertia about x
    expected_mass[4][4] = I_y;  // Rotational inertia about y
    expected_mass[5][5] = I_p;  // Polar rotational inertia

    ExpectMatrixEqual(beam_section.M_star, expected_mass);
}

}  // namespace kynema_fmb::beams::tests
