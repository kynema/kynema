#include <numbers>

#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

/**
 * @brief Static bending benchmark test for a twisted composite cantilever beam
 *
 * @details This problem demonstrates the ability of Kynema to analyze a twisted composite beam
 * with geometric twist along its length. This is an important use case for wind turbine blade
 * analysis, as turbine blades typically have significant twist from root to tip to optimize
 * aerodynamic performance.
 *
 * In this benchmark problem, we examine the static tip deflection of a straight, but twisted
 * beam. The beam has 90 degrees of twist from root to tip and is subjected to a transverse
 * tip load of 4000 kN applied in the negative z-direction. The beam is 10 units in length
 * and has uniform cross-sectional properties with isotropic material behavior (no bend-twist
 * coupling in the material, but geometric coupling due to the twist).
 *
 * The benchmark solution is a highly refined solid-element model solved in ANSYS. The problem
 * was solved with a single Legendre Spectral Finite Element (LSFE) with N-point Gauss-Legendre
 * quadrature, where N is the number of nodes. The results show the expected spectral convergence
 * to a solution that is in good agreement with the ANSYS solid-element solution. The results are
 * also identical, within machine precision, to those produced with BeamDyn.
 *
 * @note Details of this benchmark problem are described in:
 *       Q. Wang, M.A. Sprague, J. Jonkman, N. Johnson, B. Jonkman (2017).
 *       "BeamDyn: a high-fidelity wind turbine blade solver in the FAST modular framework."
 *       Wind Energy 20: 1439-1462.
 *
 * @see Kynema documentation for the full benchmark results:
 *      https://kynema.github.io/kynema/testing/twisted.html
 */
TEST(VerificationTest, CompositeBeamBending_TwistedBeam) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    // Static analysis with tight convergence tolerances
    builder.Solution()
        .EnableStaticSolve()                // Static equilibrium problem
        .SetTimeStep(1.)                    // Step size (irrelevant for static)
        .SetDampingFactor(1.)               // No numerical damping (ρ_∞ = 1)
        .SetMaximumNonlinearIterations(15)  // Max Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-11)   // Absolute tolerance
        .SetRelativeErrorTolerance(1e-9);   // Relative tolerance

    if (write_output) {
        builder.Outputs().SetOutputFilePath("StaticVerificationTest.CompositeBeamBending_TwistedBeam"
        );
    }

    //----------------------------------
    // beam element with twist
    //----------------------------------
    const int num_nodes{15};  // number of nodes = n
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // 15-node LSFE for high accuracy
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature for integration
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::Segmented)
        .PrescribedRootMotion(true);  // Root node is fixed (clamped BC)

    // Helper function to convert degrees to radians
    const auto degree_to_radians = [](const double degree) {
        return degree * std::numbers::pi / 180.;
    };

    // Define twist distribution: 0 degrees at root, 90 degrees at tip (linear variation in between)
    builder.Blade()
        .AddRefAxisTwist(0., degree_to_radians(0.))    // s = 0: twist = 0 degrees
        .AddRefAxisTwist(1., degree_to_radians(90.));  // s = 1: twist = 90 degrees

    // Reference axis geometry: straight beam along x-axis from (0,0,0) to (10,0,0)
    const std::vector<double> kp_s{0., 1.};
    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.},  // x = 10*s, y = 0, z = 0 (straight centerline)
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

    // Sectional stiffness matrix (6x6) - transversely isotropic (diagonal, no material coupling)
    constexpr auto stiffness_matrix = std::array{
        std::array{2.50e10, 0., 0., 0., 0., 0.},       //
        std::array{0., 8.2604167e9, 0., 0., 0., 0.},   //
        std::array{0., 0., 8.2604167e9, 0., 0., 0.},   //
        std::array{0., 0., 0., 1.41872656e8, 0., 0.},  //
        std::array{0., 0., 0., 0., 5.20833e8, 0.},     //
        std::array{0., 0., 0., 0., 0., 1.30208333e8},  //
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
    // Point force P_z = -4000 kN
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[2] = -4000.e3;

    // Static step
    const auto converged = interface.Step();

    // Verify convergence
    ASSERT_EQ(converged, true);

    //-------------------------------------------
    // verify tip displacements
    //-------------------------------------------
    EXPECT_NEAR(
        tip_node.displacement[0], -1.14152599953942, 1e-12
    );  // Equivalent BeamDyn soln: -1.1415259995394
        // Ansys solid element solution: -1.134192
    EXPECT_NEAR(
        tip_node.displacement[1], -1.71805334646388, 1e-12
    );  // Equivalent BeamDyn soln: 1.7180533464639
        // Ansys solid element solution: -1.714467
    EXPECT_NEAR(
        tip_node.displacement[2], -3.59321324026952, 1e-12
    );  // Equivalent BeamDyn soln: -3.5932132402695
        // Ansys solid element solution: -3.58423
}

}  // namespace kynema::tests
