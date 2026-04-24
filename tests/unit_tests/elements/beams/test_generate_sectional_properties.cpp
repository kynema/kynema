#include <array>
#include <cstddef>
#include <ranges>

#include <gtest/gtest.h>

#include "elements/beams/generate_sectional_properties.hpp"

namespace kynema_fmb::beams::tests {

class GenerateSectionalMatricesTest : public ::testing::Test {
protected:
    static constexpr double kTolerance = 1e-12;

    struct TestParameters {
        // Stiffness properties
        double EA{1000.};
        double EI_x{200.};
        double EI_y{300.};
        double GKt{150.};
        double GA{400.};
        double kxs{1.};
        double kys{1.};

        // Geometric offsets and angles
        double x_C{0.};
        double y_C{0.};
        double theta_p{0.};
        double x_S{0.};
        double y_S{0.};
        double theta_s{0.};

        // Mass properties
        double m{10.};
        double I_x{5.};
        double I_y{8.};
        double I_p{13.};  // I_x + I_y
        double x_G{0.};
        double y_G{0.};
        double theta_i{0.};
    };

    TestParameters tp_;

    [[nodiscard]] std::array<std::array<double, 6>, 6> GenerateStiffnessMatrixWithParams() const {
        return GenerateStiffnessMatrix(
            tp_.EA, tp_.EI_x, tp_.EI_y, tp_.GKt, tp_.GA, tp_.kxs, tp_.kys, tp_.x_C, tp_.y_C,
            tp_.theta_p, tp_.x_S, tp_.y_S, tp_.theta_s
        );
    }

    [[nodiscard]] std::array<std::array<double, 6>, 6> GenerateMassMatrixWithParams() const {
        return GenerateMassMatrix(tp_.m, tp_.I_x, tp_.I_y, tp_.I_p, tp_.x_G, tp_.y_G, tp_.theta_i);
    }

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

TEST_F(GenerateSectionalMatricesTest, StiffnessMatrix_Uncoupled) {
    // Test with no coupling (all offsets zero, principal angles zero)
    auto result = GenerateStiffnessMatrixWithParams();

    auto expected = std::array<std::array<double, 6>, 6>{};
    expected[0][0] = tp_.GA;    // Shear stiffness in x
    expected[1][1] = tp_.GA;    // Shear stiffness in y
    expected[2][2] = tp_.EA;    // Axial stiffness
    expected[3][3] = tp_.EI_x;  // Bending stiffness about x
    expected[4][4] = tp_.EI_y;  // Bending stiffness about y
    expected[5][5] = tp_.GKt;   // Torsional stiffness

    ExpectMatrixEqual(result, expected);
}

TEST_F(GenerateSectionalMatricesTest, StiffnessMatrix_CentroidOffsetCoupling) {
    // Test axial-bending coupling due to centroid offset
    tp_.x_C = 0.1;  // Centroid offset in x
    tp_.y_C = 0.2;  // Centroid offset in y

    auto result = GenerateStiffnessMatrixWithParams();

    // Expected stiffness matrix
    auto expected = std::array<std::array<double, 6>, 6>{};

    // Shear stiffness
    expected[0][0] = tp_.GA;
    expected[1][1] = tp_.GA;

    // Axial stiffness
    expected[2][2] = tp_.EA;

    // Torsional stiffness
    expected[5][5] = tp_.GKt;

    // Axial-bending coupling
    expected[2][3] = tp_.EA * tp_.y_C;
    expected[2][4] = -tp_.EA * tp_.x_C;
    expected[3][2] = tp_.EA * tp_.y_C;
    expected[4][2] = -tp_.EA * tp_.x_C;

    // Bending stiffness with axial coupling
    expected[3][3] = tp_.EI_x + tp_.EA * tp_.y_C * tp_.y_C;
    expected[4][4] = tp_.EI_y + tp_.EA * tp_.x_C * tp_.x_C;
    expected[3][4] = -tp_.EA * tp_.x_C * tp_.y_C;
    expected[4][3] = -tp_.EA * tp_.x_C * tp_.y_C;

    ExpectMatrixEqual(result, expected);
}

TEST_F(GenerateSectionalMatricesTest, MassMatrix_Uncoupled) {
    // Test basic mass matrix with no coupling (CG at origin, no principal axis rotation)
    auto result = GenerateMassMatrixWithParams();

    auto expected = std::array<std::array<double, 6>, 6>{};

    // Translational mass
    expected[0][0] = tp_.m;
    expected[1][1] = tp_.m;
    expected[2][2] = tp_.m;

    // Rotational inertia (no rotation, so principal = reference frame)
    expected[3][3] = tp_.I_x;
    expected[4][4] = tp_.I_y;
    expected[5][5] = tp_.I_p;

    ExpectMatrixEqual(result, expected);
}

TEST_F(GenerateSectionalMatricesTest, MassMatrix_CGOffsetCoupling) {
    // Test translational-rotational coupling due to CG offset
    tp_.x_G = 0.2;  // CG offset in x
    tp_.y_G = 0.3;  // CG offset in y

    auto result = GenerateMassMatrixWithParams();

    auto expected = std::array<std::array<double, 6>, 6>{};

    // Translational mass
    expected[0][0] = tp_.m;
    expected[1][1] = tp_.m;
    expected[2][2] = tp_.m;

    // Translational-rotational coupling due to CG offset
    expected[0][5] = -tp_.m * tp_.y_G;  // F_x due to α_z
    expected[1][5] = tp_.m * tp_.x_G;   // F_y due to α_z
    expected[2][3] = tp_.m * tp_.y_G;   // F_z due to α_x
    expected[2][4] = -tp_.m * tp_.x_G;  // F_z due to α_y

    // Rotational-translational coupling due to CG offset (symmetric)
    expected[3][2] = tp_.m * tp_.y_G;   // M_x due to a_z
    expected[4][2] = -tp_.m * tp_.x_G;  // M_y due to a_z
    expected[5][0] = -tp_.m * tp_.y_G;  // M_z due to a_x
    expected[5][1] = tp_.m * tp_.x_G;   // M_z due to a_y

    // Rotational inertia with CG coupling
    expected[3][3] = tp_.I_x + tp_.m * tp_.y_G * tp_.y_G;
    expected[4][4] = tp_.I_y + tp_.m * tp_.x_G * tp_.x_G;
    expected[3][4] = -tp_.m * tp_.x_G * tp_.y_G;
    expected[4][3] = -tp_.m * tp_.x_G * tp_.y_G;

    // Polar inertia with CG coupling
    expected[5][5] = tp_.I_p + tp_.m * (tp_.x_G * tp_.x_G + tp_.y_G * tp_.y_G);

    ExpectMatrixEqual(result, expected);
}

}  // namespace kynema_fmb::beams::tests
