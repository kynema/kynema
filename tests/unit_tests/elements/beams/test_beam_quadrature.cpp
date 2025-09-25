#include <ranges>
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

//--------------------------------------------------------------------------
// Gauss-Legendre Quadrature
//--------------------------------------------------------------------------

// Expected results were sourced from:
// https://pomax.github.io/bezierinfo/legendre-gauss.html
// and internal mathematica notebook (GaussLegendreQuadrature.nb)

TEST(BeamQuadratureTest, GaussLegendreQuadrature_1point) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(1);
    constexpr auto expected_quadrature = std::array<std::array<double, 2>, 1>{
        std::array{0., 2.}  // pt_1: 0, wt_1: 2
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_2points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(2);
    const auto expected_quadrature = std::array<std::array<double, 2>, 2>{
        std::array{-1. / std::numbers::sqrt3, 1.},  // pt_1: -1 / sqrt(3), wt_1: 1
        {1. / std::numbers::sqrt3, 1.}              // pt_2: 1 / sqrt(3), wt_2: 1
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_3points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(3);
    const auto expected_quadrature = std::array<std::array<double, 2>, 3>{
        std::array{-std::sqrt(3. / 5.), 5. / 9.},  // pt_1: -sqrt(3 / 5), wt_1: 5 / 9
        {0., 8. / 9.},                             // pt_2: 0, wt_2: 8 / 9
        {std::sqrt(3. / 5.), 5. / 9.}              // pt_3: sqrt(3 / 5), wt_3: 5 / 9
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_5points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(5);
    const auto expected_quadrature = std::array<std::array<double, 2>, 5>{
        std::array{
            -std::sqrt(5. + 2. * std::sqrt(10. / 7.)) / 3., (322. - 13. * std::sqrt(70.)) / 900.
        },
        {-std::sqrt(5. - 2. * std::sqrt(10. / 7.)) / 3., (322. + 13. * std::sqrt(70.)) / 900.},
        {0., 128. / 225.},
        {std::sqrt(5. - 2. * std::sqrt(10. / 7.)) / 3., (322. + 13. * std::sqrt(70.)) / 900.},
        {std::sqrt(5. + 2. * std::sqrt(10. / 7.)) / 3., (322. - 13. * std::sqrt(70.)) / 900.}
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_7points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(7);
    const auto expected_quadrature = std::array<std::array<double, 2>, 7>{
        std::array{-0.9491079123427585, 0.1294849661688697},
        {-0.7415311855993945, 0.2797053914892766},
        {-0.4058451513773972, 0.3818300505051189},
        {0.0, 0.4179591836734694},
        {0.4058451513773972, 0.3818300505051189},
        {0.7415311855993945, 0.2797053914892766},
        {0.9491079123427585, 0.1294849661688697}
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_9points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(9);
    const auto expected_quadrature = std::array<std::array<double, 2>, 9>{
        std::array{-0.9681602395076261, 0.0812743883615744},
        {-0.8360311073266358, 0.1806481606948574},
        {-0.6133714327005904, 0.2606106964029354},
        {-0.3242534234038089, 0.3123470770400029},
        {0.0, 0.3302393550012598},
        {0.3242534234038089, 0.3123470770400029},
        {0.6133714327005904, 0.2606106964029354},
        {0.8360311073266358, 0.1806481606948574},
        {0.9681602395076261, 0.0812743883615744}
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_10points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(10);
    // NOLINTBEGIN(modernize-use-std-numbers)
    const auto expected_quadrature = std::array<std::array<double, 2>, 10>{
        std::array{-0.973906528517172, 0.06667134430868793},
        {-0.8650633666889761, 0.1494513491504942},
        {-0.6794095682990245, 0.2190863625159827},
        {-0.4333953941292472, 0.2692667193099965},
        {-0.1488743389816312, 0.2955242247147527},
        {0.1488743389816312, 0.2955242247147527},
        {0.4333953941292472, 0.2692667193099965},
        {0.6794095682990245, 0.2190863625159827},
        {0.8650633666889761, 0.1494513491504942},
        {0.973906528517172, 0.06667134430868793}
    };
    // NOLINTEND(modernize-use-std-numbers)

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-12);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_11points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(11);
    const auto expected_quadrature = std::array<std::array<double, 2>, 11>{
        std::array{-0.9782286581460570, 0.0556685671161737},
        {-0.8870625997680953, 0.1255803694649046},
        {-0.7301520055740494, 0.1862902109277343},
        {-0.5190961292068118, 0.2331937645919905},
        {-0.2695431559523449, 0.2628045445102467},
        {0.0, 0.2729250867779006},
        {0.2695431559523449, 0.2628045445102467},
        {0.5190961292068118, 0.2331937645919905},
        {0.7301520055740494, 0.1862902109277343},
        {0.8870625997680953, 0.1255803694649046},
        {0.9782286581460570, 0.0556685671161737}
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_13points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(13);
    const auto expected_quadrature = std::array<std::array<double, 2>, 13>{
        std::array{-0.9841830547185881, 0.0404840047653159},
        {-0.9175983992229779, 0.0921214998377284},
        {-0.8015780907333099, 0.1388735102197872},
        {-0.6423493394403402, 0.1781459807619457},
        {-0.4484927510364469, 0.2078160475368885},
        {-0.2304583159551348, 0.2262831802628972},
        {0.0, 0.2325515532308739},
        {0.2304583159551348, 0.2262831802628972},
        {0.4484927510364469, 0.2078160475368885},
        {0.6423493394403402, 0.1781459807619457},
        {0.8015780907333099, 0.1388735102197872},
        {0.9175983992229779, 0.0921214998377284},
        {0.9841830547185881, 0.0404840047653159}
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

TEST(BeamQuadratureTest, GaussLegendreQuadrature_15points) {
    const auto calculated_quadrature = CreateGaussLegendreQuadrature(15);
    const auto expected_quadrature = std::array<std::array<double, 2>, 15>{
        std::array{-0.9879925180204854, 0.0307532419961173},
        {-0.9372733924007059, 0.0703660474881081},
        {-0.8482065834104272, 0.1071592204671719},
        {-0.7244177313601700, 0.1395706779261543},
        {-0.5709721726085388, 0.1662692058169939},
        {-0.3941513470775634, 0.1861610000155622},
        {-0.2011940939974345, 0.1984314853271116},
        {0.0, 0.2025782419255613},
        {0.2011940939974345, 0.1984314853271116},
        {0.3941513470775634, 0.1861610000155622},
        {0.5709721726085388, 0.1662692058169939},
        {0.7244177313601700, 0.1395706779261543},
        {0.8482065834104272, 0.1071592204671719},
        {0.9372733924007059, 0.0703660474881081},
        {0.9879925180204854, 0.0307532419961173}
    };

    for (size_t i = 0; i < expected_quadrature.size(); ++i) {
        for (size_t j = 0; j < expected_quadrature[i].size(); ++j) {
            EXPECT_NEAR(calculated_quadrature[i][j], expected_quadrature[i][j], 1e-14);
        }
    }
}

}  // namespace kynema::beams::tests
