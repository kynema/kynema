
#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"
#include "step/step.hpp"

namespace kynema::tests {

TEST(StaticVerificationTest, IsotropicBeamRollup) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    // Static analysis with tight convergence tolerances for benchmark accuracy
    builder.Solution()
        .EnableStaticSolve()
        .SetTimeStep(1.)
        .SetDampingFactor(1.)
        .SetMaximumNonlinearIterations(15)
        .SetAbsoluteErrorTolerance(1e-11)
        .SetRelativeErrorTolerance(1e-9);
    if (write_output) {
        builder.Solution().SetOutputFile("StaticVerificationTest.IsotropicBeamRollup");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{15};                      // number of nodes = n
    const int element_order{num_nodes - 1};       // 15-node LSFE for high accuracy
    const int section_refinement{num_nodes - 1};  // n-pt Gauss-Legendre quadrature
    builder.Blade()
        .SetElementOrder(element_order)
        .SetSectionRefinement(section_refinement)
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
    // Apply moment to create circular rollup
    //-------------------------------------------
    // For a beam to roll into a complete circle:
    // Curvature κ = 2π/L, Moment M = EI*κ = EI * 2π/L
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
    EXPECT_NEAR(tip_node.displacement[0], -10.0000000000645, 1e-12);  // Exact soln: -10.
    EXPECT_NEAR(tip_node.displacement[2], 0., 1e-12);                 // Exact soln: 0.
}

}  // namespace kynema::tests
