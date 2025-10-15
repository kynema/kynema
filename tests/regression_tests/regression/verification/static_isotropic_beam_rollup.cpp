
#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

/**
 * @brief Static pure bending benchmark test for an isotropic cantilever beam rolling into a circle
 *
 * @details This problem demonstrates Kynema's capability to analyze an isotropic beam with no
 * initial curvature and with highly nonlinear deflections. In fact, the displacement is so large
 * that the beam bends into a complete circular shape (2π radians of rotation over the 10-unit
 * length), providing a rigorous test of the geometrically exact beam formulation.
 *
 * The exact analytical solution for the tip displacement in the direction of the beam length
 * is -10.0, representing a complete circular rollup where the tip returns to the same
 * x-coordinate as the root.
 *
 * @note Details of this benchmark problem are described in:
 *       J. C. Simo and L. Vu-Quoc (1986). "A three-dimensional finite-strain rod model. Part II."
 *       Computer Methods in Applied Mechanics and Engineering, 58:79–116.
 *
 * @see Kynema documentation for the full benchmark results (we are using the 15-node, 15 quadrature
 *      points with Gauss-Legendre quadrature case here):
 *      https://kynema.github.io/kynema/testing/rollup.html
 */
TEST(StaticVerificationTest, IsotropicBeamRollup) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    // Static analysis with tight convergence tolerances for benchmark accuracy
    builder.Solution()
        .EnableStaticSolve()   // Static analysis
        .SetTimeStep(1.)       // Step size (irrelevant for static)
        .SetDampingFactor(1.)  // No numerical damping (ρ_∞ = 1, irrelevant for static)
        .SetMaximumNonlinearIterations(15)  // Max Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-11)
        .SetRelativeErrorTolerance(1e-9);
    if (write_output) {
        builder.Solution().SetOutputFile("StaticVerificationTest.IsotropicBeamRollup");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{15};  // number of nodes = n
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // 15-node LSFE for high accuracy
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature for integration
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::WholeBeam)
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
    const std::vector<double> section_s{0., 1.};
    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            // mass matrix -- just a placeholder for static analysis
            std::array{
                std::array{1., 0., 0., 0., 0., 0.},
                std::array{0., 1., 0., 0., 0., 0.},
                std::array{0., 0., 1., 0., 0., 0.},
                std::array{0., 0., 0., 1., 0., 0.},
                std::array{0., 0., 0., 0., 1., 0.},
                std::array{0., 0., 0., 0., 0., 1.},
            },
            // stiffness matrix -- isotropic beam
            std::array{
                std::array{1770.e3, 0., 0., 0., 0., 0.},
                std::array{0., 1770.e3, 0., 0., 0., 0.},
                std::array{0., 0., 1770.e3, 0., 0., 0.},
                std::array{0., 0., 0., 8.16e3, 0., 0.},
                std::array{0., 0., 0., 0., 86.9e3, 0.},
                std::array{0., 0., 0., 0., 0., 215.e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    //-------------------------------------------
    // apply moment to create circular rollup
    //-------------------------------------------
    // For a beam to roll into a complete circle -> curvature k = 2π/L, moment = EI * k = EI * 2π/L
    const auto moment = 2. * M_PI * 86.9e3 / 10.;

    // Apply moment to tip node about y axis (negative for rollup)
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[4] = -moment;

    // Static step
    const auto converged = interface.Step();

    // Check convergence
    ASSERT_EQ(converged, true);

    //----------------------------------
    // verify tip displacements
    //----------------------------------
    EXPECT_NEAR(tip_node.displacement[0], -10.0000000000645, 1e-12);  // Exact analytical soln: -10.
    EXPECT_NEAR(tip_node.displacement[2], 0., 1e-12);                 // Exact analytical soln: 0.
}

}  // namespace kynema::tests
