#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

TEST(StaticVerificationTest, CompositeBeamBending) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    // Static analysis with tight convergence tolerances
    builder.Solution()
        .EnableStaticSolve()   // Static analysis
        .SetTimeStep(1.)       // Step size (irrelevant for static)
        .SetDampingFactor(1.)  // No numerical damping (ρ_∞ = 1, irrerelevant for static)
        .SetMaximumNonlinearIterations(15)  // Max Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-11)
        .SetRelativeErrorTolerance(1e-9);

    if (write_output) {
        builder.Solution().SetOutputFile("StaticVerificationTest.CompositeBeamBending");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{15};  // 15-node LSFE
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // 15-node LSFE for high accuracy
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::WholeBeam)
        .PrescribedRootMotion(true);  // Root node is fixed (clamped BC)

    // No twist along beam reference axis (straight beam)
    builder.Blade()
        .AddRefAxisTwist(0., 0.)   // s = 0: twist = 0
        .AddRefAxisTwist(1., 0.);  // s = 1: twist = 0

    // Reference axis geometry: straight beam along x-axis from (0,0,0) to (10,0,0)
    const std::vector<double> kp_s{0., 1.};
    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.},  // x = 10*s, y = 0, z = 0 (straight beam)
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    //----------------------------------
    // beam cross-section properties
    //----------------------------------

    // Mass matrix
    constexpr auto mass_matrix = std::array{
        std::array{8.538e-2, 0., 0., 0., 0., 0.},   std::array{0., 8.538e-2, 0., 0., 0., 0.},
        std::array{0., 0., 8.538e-2, 0., 0., 0.},   std::array{0., 0., 0., 1.4433e-2, 0., 0.},
        std::array{0., 0., 0., 0., 0.40972e-2, 0.}, std::array{0., 0., 0., 0., 0., 1.0336e-2},
    };

    // Stiffness matrix
    constexpr auto stiffness_matrix = std::array{
        std::array{1368.17e3, 0., 0., 0., 0., 0.},
        std::array{0., 88.56e3, 0., 0., 0., 0.},
        std::array{0., 0., 38.78e3, 0., 0., 0.},
        std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
        std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
        std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
    };

    // Apply uniform properties along entire beam length
    const std::vector<double> section_s{0., 1.};
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s, mass_matrix, stiffness_matrix, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    //-------------------------------------------
    // Apply transverse tip load
    //-------------------------------------------
    // Point force P_z = 150 lbs
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[2] = 150.;

    // Static step
    const auto converged = interface.Step();

    // Verify convergence
    ASSERT_EQ(converged, true);

    //-------------------------------------------
    // Verify tip displacements
    //-------------------------------------------
    EXPECT_NEAR(
        tip_node.displacement[0], -9.02726627566299E-02, 1e-12
    );  // Beamdyn soln: -0.09064180058
    EXPECT_NEAR(
        tip_node.displacement[1], -6.47488486259036E-02, 1e-12
    );  // Beamdyn soln: -0.0648265579
    EXPECT_NEAR(tip_node.displacement[2], 1.22973648292371E+00, 1e-12);  // Beamdyn soln: 1.229985112
}

}  // namespace kynema::tests
