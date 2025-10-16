#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

/**
 * @brief Static bending benchmark test for a curved composite cantilever beam
 *
 * @details This problem demonstrates the ability of Kynema to analyze beams with initial
 * curvature in their reference configuration. This is an important capability for wind
 * turbine blade analysis, as many modern blade designs incorporate prebend (out-of-plane
 * curvature) and sweep (in-plane curvature) to optimize structural and aerodynamic
 * performance.
 *
 * In this benchmark problem, we examine the tip deflection of a tip-loaded curved beam.
 * The beam lies in the x-y plane (initially) with a tip force of 600 lbs acting in the
 * positive z direction. The beam is curved on a radius of 100 inches in the positive x
 * direction through a 45-degree arc, creating significant initial curvature that couples
 * with the applied loads to produce complex 3D deformations.
 *
 * @note Details of this benchmark problem are described in:
 *       K.-J. Bathe and S. Bolourchi (1979).
 *       "Large displacement analysis of three-dimensional beam structures."
 *       International Journal for Numerical Methods in Engineering 14: 961-986.
 *
 * @see Kynema documentation for the full benchmark results and convergence study:
 *      https://kynema.github.io/kynema/testing/curved.html
 */
TEST(StaticVerificationTest, CompositeBeamBending_CurvedBeam) {
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
        .SetAbsoluteErrorTolerance(1e-11)   // absolute tolerance
        .SetRelativeErrorTolerance(1e-9);   // relative tolerance

    if (write_output) {
        builder.Solution().SetOutputFile("StaticVerificationTest.CompositeBeamBending_CurvedBeam");
    }

    //----------------------------------
    // beam element w/ curved geometry
    //----------------------------------
    const int num_nodes{15};  // number of nodes = n
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // 15-node LSFE for high accuracy
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature for integration
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::Segmented)
        .PrescribedRootMotion(true);  // Root node is fixed (clamped BC)

    // No twist along beam reference axis (planar curvature only)
    builder.Blade()
        .AddRefAxisTwist(0., 0.)   // s = 0: twist = 0°
        .AddRefAxisTwist(1., 0.);  // s = 1: twist = 0°

    // Define curved reference axis geometry (circular arc in x-y plane)
    // Arc spans approximately 45 degrees of a circle with radius ~100
    const std::vector<std::tuple<double, double, double>> arc_points{
        {0., 0., 0.},                     // Root at origin
        {3.82683e+01, -7.61205e+00, 0.},  // Mid-arc point
        {7.07107e+01, -2.92893e+01, 0.}   // Tip (approximately at 45° angle)
    };
    for (const auto& s : arc_points) {
        builder.Blade().AddRefAxisPoint(
            std::get<0>(s), {std::get<0>(s), std::get<1>(s), std::get<2>(s)},
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

    // Sectional stiffness matrix (6x6) - isotropic (diagonal, no coupling)
    constexpr auto stiffness_matrix = std::array{
        std::array{100.e5, 0., 0., 0., 0., 0.},     //
        std::array{0., 41.6667e5, 0., 0., 0., 0.},  //
        std::array{0., 0., 41.6667e5, 0., 0., 0.},  //
        std::array{0., 0., 0., 8.333e5, 0., 0.},    //
        std::array{0., 0., 0., 0., 8.333e5, 0.},    //
        std::array{0., 0., 0., 0., 0., 8.333e5},    //
    };

    // Apply uniform properties along entire beam arc length
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
    // Point force P_z = 600 in positive z-direction
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[2] = 600.;

    // Static step
    const auto converged = interface.Step();

    // Verify convergence
    ASSERT_EQ(converged, true);

    //-------------------------------------------
    // verify tip displacements
    //-------------------------------------------
    EXPECT_NEAR(
        tip_node.displacement[0], -23.6459716009047, 1e-12
    );  // Equivalent BeamDyn soln: -23.645971600905
        // Published solution: -23.5
    EXPECT_NEAR(
        tip_node.displacement[1], 13.7208006040058, 1e-12
    );  // Equivalent BeamDyn soln: -13.720800604006
        // Published solution: 13.4
    EXPECT_NEAR(
        tip_node.displacement[2], 53.5776450626987, 1e-12
    );  // Equivalent BeamDyn soln: 53.577645062699
        // Published solution: 53.4
}

}  // namespace kynema::tests
