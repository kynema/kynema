#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

/**
 * @brief Static bending benchmark test for an anisotropic, composite cantilever beam with coupling
 * effects
 *
 * @details This problem demonstrates the ability of Kynema to analyze a composite beam with
 * bend-bend and bend-twist coupling. This is an important use case for Kynema to be able to model
 * since virtually all wind turbine blades are constructed with such coupling effects, where
 * off-diagonal terms in the 6x6 sectional stiffness matrix create interactions between bending
 * and torsional modes.
 *
 * In this benchmark problem, we examine the tip deflection of a composite beam under a tip load
 * of 150 lbs applied in the z-direction. The beam is 10 units in length and has uniform
 * cross-sectional properties consisting of anistropic material properties. The problem was solved
 * with a single Legendre Spectral Finite Element (LSFE) with N-point Gauss-Legendre quadrature where
 * N is the number of nodes.
 *
 * @note Details of this benchmark problem are described in:
 *       Q. Wang, M.A. Sprague, J. Jonkman, N. Johnson, B. Jonkman (2017).
 *       "BeamDyn: a high-fidelity wind turbine blade solver in the FAST modular framework."
 *       Wind Energy 20: 1439-1462.
 *
 * @see Kynema documentation for the full benchmark results (we are using the 15-node, 15 quadrature
 *      points with Gauss-Legendre quadrature case here):
 *      https://kynema.github.io/kynema/testing/composite.html
 */
TEST(VerificationTest, Static_CompositeBeamBending) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    // Static analysis with tight convergence tolerances for benchmark accuracy
    builder.Solution()
        .EnableStaticSolve()   // Static analysis
        .SetTimeStep(1.)       // Step size (irrelevant for static)
        .SetDampingFactor(1.)  // No numerical damping (ρ_∞ = 1, irrerelevant for static)
        .SetMaximumNonlinearIterations(15)  // Max Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-11)
        .SetRelativeErrorTolerance(1e-9);

    if (write_output) {
        builder.Outputs().SetOutputFilePath("VerificationTest.Static_CompositeBeamBending");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{15};  // number of nodes = n
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // 15-node LSFE for high accuracy
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature for integration
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::Segmented)
        .PrescribedRootMotion(true);  // Root node is fixed (clamped BC)

    // No twist along beam reference axis
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

    // Sectional mass matrix (6x6) - placeholder for static analysis
    constexpr auto mass_matrix = std::array{
        std::array{1., 0., 0., 0., 0., 0.},  //
        std::array{0., 1., 0., 0., 0., 0.},  //
        std::array{0., 0., 1., 0., 0., 0.},  //
        std::array{0., 0., 0., 1., 0., 0.},  //
        std::array{0., 0., 0., 0., 1., 0.},  //
        std::array{0., 0., 0., 0., 0., 1.},  //
    };

    // sectional stiffness matrix (6x6)
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
    // apply transverse tip load
    //-------------------------------------------
    // Point force P_z = 150 lbs
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[2] = 150.;

    // Static step
    const auto converged = interface.Step();

    // Verify convergence
    ASSERT_EQ(converged, true);

    //-------------------------------------------
    // verify tip displacements
    //-------------------------------------------
    EXPECT_NEAR(
        tip_node.displacement[0], -9.02726627566299E-02, 1e-12
    );  // Equivalent BeamDyn soln: -0.090272662756631
    EXPECT_NEAR(
        tip_node.displacement[1], -6.47488486259036E-02, 1e-12
    );  // Equivalent BeamDyn soln: -0.064748848625905
    EXPECT_NEAR(
        tip_node.displacement[2], 1.22973648292371E+00, 1e-12
    );  // Equivalent BeamDyn soln: 1.2297364829237
}

}  // namespace kynema::tests
