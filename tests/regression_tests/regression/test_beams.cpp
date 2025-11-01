#include <array>
#include <initializer_list>
#include <vector>

#include <gtest/gtest.h>

#include "model/model.hpp"
#include "test_utilities.hpp"

namespace {
inline auto SetUpBeams() {
    auto model = kynema::Model();

    model.SetGravity(0., 0., 9.81);

    // Add nodes and get list of node IDs
    const std::vector<size_t> beam_node_ids({
        model.AddNode()
            .SetElemLocation(0.)
            .SetPosition(
                0., 0., 0, 0.9778215200524469, -0.01733607539094763, -0.09001900002195001,
                -0.18831121859148398
            )
            .Build(),
        model.AddNode()
            .SetElemLocation(0.1726731646460114)
            .SetPosition(
                0.863365823230057, -0.2558982639254171, 0.11304112106827427, 0.9950113028068008,
                -0.002883848832932071, -0.030192109815745303, -0.09504013471947484
            )
            .SetDisplacement(
                0.002981602178886856, -0.00246675949494302, 0.003084570715675624, 0.9999627302042724,
                0.008633550973807708, 0., 0.
            )
            .SetVelocity(
                0.01726731646460114, -0.014285714285714285, 0.003084570715675624,
                0.01726731646460114, -0.014285714285714285, 0.003084570715675624
            )
            .SetAcceleration(
                0.01726731646460114, -0.011304112106827427, 0.00606617289456248, 0.01726731646460114,
                -0.014285714285714285, -0.014285714285714285
            )
            .Build(),
        model.AddNode()
            .SetElemLocation(0.5)
            .SetPosition(
                2.5, -0.25, 0, 0.9904718430204884, -0.009526411091536478, 0.09620741150793366,
                0.09807604012323785
            )
            .SetDisplacement(
                0.025, -0.0125, 0.027500000000000004, 0.9996875162757026, 0.02499739591471221, 0, 0
            )
            .SetVelocity(0.05, -0.025, 0.027500000000000004, 0.05, -0.025, 0.027500000000000004)
            .SetAcceleration(0.05, 0, 0.052500000000000005, 0.05, -0.025, -0.025)
            .Build(),
        model.AddNode()
            .SetElemLocation(0.82732683535398865)
            .SetPosition(
                4.1366341767699435, 0.39875540678256005, -0.5416125496397031, 0.9472312341234699,
                -0.04969214162931507, 0.18127630174800594, 0.25965858850765167
            )
            .SetDisplacement(
                0.06844696924968459, -0.011818954790771264, 0.07977257214146725, 0.9991445348823055,
                0.04135454527402512, 0., 0.
            )
            .SetVelocity(
                0.08273268353539887, -0.01428571428571428, 0.07977257214146725, 0.08273268353539887,
                -0.01428571428571428, 0.07977257214146725
            )
            .SetAcceleration(
                0.08273268353539887, 0.05416125496397031, 0.14821954139115184, 0.08273268353539887,
                -0.01428571428571428, -0.01428571428571428
            )
            .Build(),
        model.AddNode()
            .SetElemLocation(1.)
            .SetPosition(
                5., 1., -1, 0.9210746582719719, -0.07193653093139739, 0.20507529985516368,
                0.32309554437664584
            )
            .SetDisplacement(0.1, 0., 0.12, 0.9987502603949663, 0.04997916927067825, 0., 0.)
            .SetVelocity(0.1, 0., 0.12, 0.1, 0., 0.12)
            .SetAcceleration(0.1, 0.1, 0.22, 0.1, 0., 0.)
            .Build(),
    });

    // Stiffness matrix for uniform composite beam section
    constexpr auto stiffness_matrix = std::array{
        std::array{1., 2., 3., 4., 5., 6.},       //
        std::array{2., 4., 6., 8., 10., 12.},     //
        std::array{3., 6., 9., 12., 15., 18.},    //
        std::array{4., 8., 12., 16., 20., 24.},   //
        std::array{5., 10., 15., 20., 25., 30.},  //
        std::array{6., 12., 18., 24., 30., 36.},  //
    };

    // Mass matrix for uniform composite beam section
    constexpr auto mass_matrix = std::array{
        std::array{2., 0., 0., 0., 0.6, -0.4},  //
        std::array{0., 2., 0., -0.6, 0., 0.2},  //
        std::array{0., 0., 2., 0.4, -0.2, 0.},  //
        std::array{0., -0.6, 0.4, 1., 2., 3.},  //
        std::array{0.6, 0., -0.2, 2., 4., 6.},  //
        std::array{-0.4, 0.2, 0., 3., 6., 9.},  //
    };

    model.AddBeamElement(
        beam_node_ids,
        std::array{
            kynema::BeamSection(0., mass_matrix, stiffness_matrix),
            kynema::BeamSection(1., mass_matrix, stiffness_matrix),
        },
        std::array{
            std::array{-0.9491079123427585, 0.1294849661688697},
            std::array{-0.7415311855993943, 0.27970539148927664},
            std::array{-0.40584515137739696, 0.3818300505051189},
            std::array{6.123233995736766e-17, 0.4179591836734694},
            std::array{0.4058451513773971, 0.3818300505051189},
            std::array{0.7415311855993945, 0.27970539148927664},
            std::array{0.9491079123427585, 0.1294849661688697},
        }
    );

    auto beams = model.CreateBeams<
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>>(
    );

    return beams;
}

}  // namespace
namespace kynema::tests {

TEST(BeamsTest, NodeInitialPositionX0) {
    const auto beams = SetUpBeams();
    expect_kokkos_view_2D_equal(
        Kokkos::subview(beams.node_x0, 0, Kokkos::ALL, Kokkos::ALL),
        {
            {0., 0., 0., 0.9778215200524469, -0.01733607539094763, -0.09001900002195001,
             -0.18831121859148398},
            {0.863365823230057, -0.2558982639254171, 0.11304112106827427, 0.9950113028068008,
             -0.002883848832932071, -0.030192109815745303, -0.09504013471947484},
            {2.5, -0.25, 0., 0.9904718430204884, -0.009526411091536478, 0.09620741150793366,
             0.09807604012323785},
            {4.1366341767699435, 0.39875540678256005, -0.5416125496397031, 0.9472312341234699,
             -0.04969214162931507, 0.18127630174800594, 0.25965858850765167},
            {5., 1., -1., 0.9210746582719719, -0.07193653093139739, 0.20507529985516368,
             0.32309554437664584},
        }
    );
}

TEST(BeamsTest, QuadraturePointMassMatrixInMaterialFrame) {
    const auto beams = SetUpBeams();
    auto Mstar = Kokkos::View<double[6][6]>("Mstar");
    Kokkos::deep_copy(Mstar, Kokkos::subview(beams.qp_Mstar, 0, 0, Kokkos::ALL, Kokkos::ALL));
    expect_kokkos_view_2D_equal(
        Mstar,
        {
            {2., 0., 0., 0., 0.6, -0.4},
            {0., 2., 0., -0.6, 0., 0.2},
            {0., 0., 2., 0.4, -0.2, 0.},
            {0., -0.6, 0.4, 1., 2., 3.},
            {0.6, 0., -0.2, 2., 4., 6.},
            {-0.4, 0.2, 0., 3., 6., 9.},

        }
    );
}

TEST(BeamsTest, QuadraturePointStiffnessMatrixInMaterialFrame) {
    const auto beams = SetUpBeams();
    auto Cstar = Kokkos::View<double[6][6]>("Cstar");
    Kokkos::deep_copy(Cstar, Kokkos::subview(beams.qp_Cstar, 0, 0, Kokkos::ALL, Kokkos::ALL));
    expect_kokkos_view_2D_equal(
        Cstar,
        {
            {1., 2., 3., 4., 5., 6.},
            {2., 4., 6., 8., 10., 12.},
            {3., 6., 9., 12., 15., 18.},
            {4., 8., 12., 16., 20., 24.},
            {5., 10., 15., 20., 25., 30.},
            {6., 12., 18., 24., 30., 36.},
        }
    );
}

}  // namespace kynema::tests
