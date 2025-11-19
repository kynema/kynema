#include <gtest/gtest.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

TEST(VerificationTest, Dynamic_ClampedCompositeBeamBending) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const double time_step{0.005};
    const double simulation_time{1.};
    builder.Solution()
        .EnableDynamicSolve()               // Dynamic analysis
        .SetTimeStep(time_step)             // Time step size
        .SetDampingFactor(0.)               // Max numerical damping (ρ_∞ = 0.)
        .SetMaximumNonlinearIterations(15)  // Max number of Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-7)    // Absolute error tolerance
        .SetRelativeErrorTolerance(1e-5);   // Relative error tolerance

    const auto write_output{false};
    if (write_output) {
        builder.Outputs().SetOutputFilePath("VerificationTest.Dynamic_ClampedCompositeBeamBending");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{5};  // number of nodes = n
    builder.Blade()
        .SetElementOrder(num_nodes - 1)       // element order = n-1
        .SetSectionRefinement(num_nodes - 1)  // n-pt Gauss-Legendre quadrature for integration
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::WholeBeam)
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
        std::array{8.538e-2, 0., 0., 0., 0., 0.},   //
        std::array{0., 8.538e-2, 0., 0., 0., 0.},   //
        std::array{0., 0., 8.538e-2, 0., 0., 0.},   //
        std::array{0., 0., 0., 1.4433e-2, 0., 0.},  //
        std::array{0., 0., 0., 0., 0.4097e-2, 0.},  //
        std::array{0., 0., 0., 0., 0., 1.0336e-2},  //
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

    //--------------------------------------------------------------------
    // apply transverse tip load and run simulation for a time period
    //--------------------------------------------------------------------
    // Point force P_z = 150 lbs applied at tip node
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];
    tip_node.loads[2] = 150.;
    const auto num_steps = static_cast<size_t>(std::round(simulation_time / time_step));

    // Reference results from BeamDyn simulation stored in: dynamic_beam_5nodes_0pt005s_beamdyn.csv
    // Comparing only the z-component (index 2) i.e. displacement along the load direction
    const std::map<double, double> beamdyn_displacements{
        {0.05, 7.6836577300452E-01},  // time = 0.05s
        {0.1, 2.0534850126406E+00},   // time = 0.10s
        {0.15, 2.2052968852838E+00},  // time = 0.15s
        {0.2, 1.1064712884566E+00},   // time = 0.20s
        {0.25, 1.5758876087082E-01},  // time = 0.25s
        {0.3, 5.5309538849698E-01},   // time = 0.30s
        {0.35, 1.8159069869944E+00},  // time = 0.35s
        {0.4, 2.2654942795857E+00},   // time = 0.40s
        {0.45, 1.3113534620467E+00},  // time = 0.45s
        {0.5, 2.7880849115281E-01},   // time = 0.50s
        {0.55, 4.4070037636959E-01},  // time = 0.55s
        {0.6, 1.5982497943061E+00},   // time = 0.60s
        {0.65, 2.2970328612154E+00},  // time = 0.65s
        {0.7, 1.5343329056122E+00},   // time = 0.70s
        {0.75, 3.4515990609276E-01},  // time = 0.75s
        {0.8, 2.7817188933310E-01},   // time = 0.80s
        {0.85, 1.3946588334800E+00},  // time = 0.85s
        {0.9, 2.2886913889084E+00},   // time = 0.90s
        {0.95, 1.8093410364286E+00},  // time = 0.95s
        {1., 5.4135602327642E-01}     // time = 1.00s
    };

    // Extract verification times from reference data
    std::vector<double> verification_times;
    verification_times.reserve(beamdyn_displacements.size());
    for (const auto& [time, _] : beamdyn_displacements) {
        verification_times.push_back(time);
    }

    // Map step numbers to verification times for exact key matching
    std::map<size_t, double> step_to_time;
    for (const auto& verify_time : verification_times) {
        const auto step = static_cast<size_t>(std::round(verify_time / time_step));
        step_to_time[step] = verify_time;
    }

    // Store tip displacements at verification time steps
    std::map<double, std::array<double, 3>> tip_displacements;
    for (auto step : std::views::iota(1U, num_steps + 1)) {
        // Take a single time step in dynamic solve
        auto converged = interface.Step();

        // Verify we reach convergence
        ASSERT_EQ(converged, true);

        // Store tip displacement at verification time steps
        if (auto it = step_to_time.find(step); it != step_to_time.end()) {
            tip_displacements[it->second] = {
                tip_node.displacement[0], tip_node.displacement[1], tip_node.displacement[2]
            };
        }
    }

    //-------------------------------------------
    // verify tip displacements
    //-------------------------------------------
    // Using relative tolerance (20%) with minimum absolute threshold for small values
    // NOTE: The tolerance is high since we don't expect Kynema's constraint-based approach
    // to match BeamDyn's Dirichlet BCs based approach
    constexpr double relative_tolerance{0.20};  // 20% relative tolerance
    constexpr double min_abs_tolerance{1e-5};   // Minimum absolute tolerance for near-zero values
    auto expect_near_with_relative_tolerance = [&](const auto& actual, const auto& expected) {
        const auto tolerance = std::max(relative_tolerance * std::abs(expected), min_abs_tolerance);
        EXPECT_NEAR(actual, expected, tolerance);
    };

    // Verify all expected times are present and compare displacements
    for (const auto& [time, ref_disp] : beamdyn_displacements) {
        const auto it = tip_displacements.find(time);
        ASSERT_NE(it, tip_displacements.end()) << "Missing tip displacement data for time: " << time;
        expect_near_with_relative_tolerance(it->second[2], ref_disp);
    }
}

}  // namespace kynema::tests
