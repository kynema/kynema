#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

TEST(DynamicVerificationTest, ClampedCompositeBeamBendingUnderTipLoad) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{true};
    const double time_step{0.01};
    builder.Solution()
        .EnableDynamicSolve()               // Dynamic analysis
        .SetTimeStep(time_step)             // Time step size
        .SetDampingFactor(0.9)              // Small numerical damping (ρ_∞ = 0.9)
        .SetMaximumNonlinearIterations(15)  // Max number of Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-7)    // Absolute error tolerance
        .SetRelativeErrorTolerance(1e-5);   // Relative error tolerance

    if (write_output) {
        builder.Outputs().SetOutputFilePath(
            "DynamicVerificationTest.ClampedCompositeBeamBendingUnderTipLoad"
        );
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{18};  // number of nodes = n
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // 18-node LSFE for high accuracy
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature for integration
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::Segmented)
        .PrescribedRootMotion(true);  // Root node is fixed (i.e. clamped BC)

    // No twist along beam reference axis
    builder.Blade()
        .AddRefAxisTwist(0., 0.)   // s = 0: twist = 0
        .AddRefAxisTwist(1., 0.);  // s = 1: twist = 0

    // Reference axis geometry: straight beam along x-axis from (0,0,0) -> (10,0,0)
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

    // Sectional mass matrix (6x6)
    constexpr auto mass_matrix = std::array{
        std::array{8.538e-2, 0., 0., 0., 0., 0.},    //
        std::array{0., 8.538e-2, 0., 0., 0., 0.},    //
        std::array{0., 0., 8.538e-2, 0., 0., 0.},    //
        std::array{0., 0., 0., 1.44332e-2, 0., 0.},  //
        std::array{0., 0., 0., 0., 0.40972e-2, 0.},  //
        std::array{0., 0., 0., 0., 0., 1.0336e-2},   //
    };

    // Sectional stiffness matrix (6x6)
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

    //------------------------------------------------------
    // apply transverse tip load and run for 10s
    //------------------------------------------------------
    // Point force P_z = 150 lbs
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[2] = 150.;

    /*
    std::cout << "Time, Tip node displacement in x direction, Tip node displacement in y direction, "
              << "Tip node displacement in z direction" << "\n";
    std::cout << 0. << ", " << std::setprecision(15) << tip_node.displacement[0] << ", "
              << tip_node.displacement[1] << ", " << tip_node.displacement[2] << "\n";
    */
    const auto num_steps = static_cast<size_t>(1. / time_step);  // 1 s at time step size = 0.01 s
    for ([[maybe_unused]] auto step : std::views::iota(1U, num_steps + 1)) {
        // Take a single time step in dynamic solve
        auto converged = interface.Step();

        // Verify we reach convergence
        ASSERT_EQ(converged, true);

        /*
        auto time = static_cast<double>(step) * time_step;
        std::cout << time << ", " << std::setprecision(15) << tip_node.displacement[0] << ", "
                  << tip_node.displacement[1] << ", " << tip_node.displacement[2] << "\n";
        */
    }

    //-------------------------------------------
    // verify tip displacements
    //-------------------------------------------

    // TODO
}

}  // namespace kynema::tests
