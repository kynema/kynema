#include <vector>

#include <gtest/gtest.h>

#include "elements/beams/beam_quadrature.hpp"

namespace kynema::beams::tests {

//--------------------------------------------------------------------------
// Trapezoidal Quadrature
//--------------------------------------------------------------------------

TEST(BeamQuadratureTest, CheckCreateTrapezoidalQuadrature_1) {
    const auto calculated_quadrature =
        CreateTrapezoidalQuadrature(std::array{0., 0.2, 0.4, 0.6, 0.8, 1.});
    constexpr auto expected_quadrature = std::array<std::array<double, 2>, 6>{
        std::array{-1., 0.2},  //
        {-0.6, 0.4},           //
        {-0.2, 0.4},           //
        {0.2, 0.4},            //
        {0.6, 0.4},            //
        {1., 0.2}              //
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, CheckCreateTrapezoidalQuadrature_2) {
    const auto calculated_quadrature =
        CreateTrapezoidalQuadrature(std::array{-5., -3., -1., 0., 3., 4., 5.});
    constexpr auto expected_quadrature = std::array<std::array<double, 2>, 7>{
        std::array{-1., 0.2},  //
        {-0.6, 0.4},           //
        {-0.2, 0.3},           //
        {0.0, 0.4},            //
        {0.6, 0.4},            //
        {0.8, 0.2},            //
        {1.0, 0.1}             //
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

}  // namespace kynema::beams::tests
