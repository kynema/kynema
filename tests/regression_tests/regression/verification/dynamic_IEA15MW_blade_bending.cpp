#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

TEST(VerificationTest, Dynamic_Iea15MwBladeBending) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const double time_step{0.001};
    const double simulation_time{1.};
    builder.Solution()
        .EnableDynamicSolve()               // Dynamic analysis
        .SetTimeStep(time_step)             // Time step size
        .SetDampingFactor(0.4)              // Medium numerical damping (ρ_∞ = 0.4)
        .SetMaximumNonlinearIterations(15)  // Max Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-5)    // Absolute error tolerance
        .SetRelativeErrorTolerance(1e-5);   // Relative error tolerance

    const auto write_output{false};
    if (write_output) {
        builder.Outputs().SetOutputFilePath("VerificationTest.Dynamic_Iea15MwBladeBending");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{11};                 // number of nodes = n
    const int element_order{num_nodes - 1};  // element order = n-1
    const int section_refinement{1};         // integrate each section w/ n-pt G-L quadrature

    builder.Blade()
        .SetElementOrder(element_order)
        .SetSectionRefinement(section_refinement)
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::Segmented)
        .PrescribedRootMotion(true);  // clamped at root

    //----------------------------------------
    // build blade w/ definition from windIO
    //----------------------------------------
    const YAML::Node wio = YAML::LoadFile("interfaces_test_files/IEA-15-240-RWT.yaml");
    const auto& wio_blade = wio["components"]["blade"];

    // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
    const auto ref_axis = wio_blade["outer_shape_bem"]["reference_axis"];
    const auto axis_grid = ref_axis["x"]["grid"].as<std::vector<double>>();
    const auto x_values = ref_axis["x"]["values"].as<std::vector<double>>();
    const auto y_values = ref_axis["y"]["values"].as<std::vector<double>>();
    const auto z_values = ref_axis["z"]["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, axis_grid.size())) {
        builder.Blade().AddRefAxisPoint(
            axis_grid[i], {x_values[i], y_values[i], z_values[i]},
            interfaces::components::ReferenceAxisOrientation::Z
        );
    }

    // Add reference axis twist
    const auto twist = wio_blade["outer_shape_bem"]["twist"];
    const auto twist_grid = twist["grid"].as<std::vector<double>>();
    const auto twist_values = twist["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, twist_grid.size())) {
        builder.Blade().AddRefAxisTwist(twist_grid[i], -twist_values[i]);
    }

    // Add blade section properties (mass and stiffness)
    const auto stiff_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["stiff_matrix"];
    const auto inertia_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["inertia_matrix"];
    const auto k_grid = stiff_matrix["grid"].as<std::vector<double>>();
    const auto m_grid = inertia_matrix["grid"].as<std::vector<double>>();
    const auto n_sections = k_grid.size();
    if (m_grid.size() != k_grid.size()) {
        throw std::runtime_error("stiffness and mass matrices not on same grid");
    }
    for (auto i : std::views::iota(0U, n_sections)) {
        if (std::abs(m_grid[i] - k_grid[i]) > 1e-8) {
            throw std::runtime_error("stiffness and mass matrices not on same grid");
        }
        const auto m = inertia_matrix["values"][i].as<std::vector<double>>();
        const auto k = stiff_matrix["values"][i].as<std::vector<double>>();
        builder.Blade().AddSection(
            m_grid[i],
            {{
                {m[0], m[1], m[2], m[3], m[4], m[5]},
                {m[1], m[6], m[7], m[8], m[9], m[10]},
                {m[2], m[7], m[11], m[12], m[13], m[14]},
                {m[3], m[8], m[12], m[15], m[16], m[17]},
                {m[4], m[9], m[13], m[16], m[18], m[19]},
                {m[5], m[10], m[14], m[17], m[19], m[20]},
            }},
            {{
                {k[0], k[1], k[2], k[3], k[4], k[5]},
                {k[1], k[6], k[7], k[8], k[9], k[10]},
                {k[2], k[7], k[11], k[12], k[13], k[14]},
                {k[3], k[8], k[12], k[15], k[16], k[17]},
                {k[4], k[9], k[13], k[16], k[18], k[19]},
                {k[5], k[10], k[14], k[17], k[19], k[20]},
            }},
            interfaces::components::ReferenceAxisOrientation::Z
        );
    }

    auto interface = builder.Build();

    //------------------------------------------------------------
    // apply load and advance simulation for a time period
    //------------------------------------------------------------
    // Apply tip load in flapwise direction i.e. -z axis, enough to cause ~10% deflection of tip node
    auto& tip_node = interface.Blade().nodes.back();
    tip_node.loads[2] = -200.e3;  // -200 kN
    const auto num_steps = static_cast<size_t>(std::round(simulation_time / time_step));

    // Reference results from BeamDyn simulation stored in:
    // IEA15blade_11nodes_rho_inf_0pt40_DT_0pt001s_F_-200kN_beamdyn.csv These values are tip
    // out-of-plane displacements (z direction) at the tip node, sampled every 0.05s from t=0.05s to
    // t=1.00s
    const std::map<double, double> beamdyn_displacements{
        {0.05, 3.2152402220957E+00},  // time = 0.05s
        {0.10, 6.2794790513560E+00},  // time = 0.10s
        {0.15, 8.4217275168128E+00},  // time = 0.15s
        {0.20, 1.0223433394254E+01},  // time = 0.20s
        {0.25, 1.1546525479042E+01},  // time = 0.25s
        {0.30, 1.2886325711343E+01},  // time = 0.30s
        {0.35, 1.3886527167187E+01},  // time = 0.35s
        {0.40, 1.4310036443370E+01},  // time = 0.40s
        {0.45, 1.6185820024709E+01},  // time = 0.45s
        {0.50, 1.7668197552712E+01},  // time = 0.50s
        {0.55, 1.6540424020922E+01},  // time = 0.55s
        {0.60, 1.5851104422733E+01},  // time = 0.60s
        {0.65, 1.5399339110826E+01},  // time = 0.65s
        {0.70, 1.6837434957033E+01},  // time = 0.70s
        {0.75, 2.0198188171276E+01},  // time = 0.75s
        {0.80, 2.2050623695122E+01},  // time = 0.80s
        {0.85, 2.4902345492894E+01},  // time = 0.85s
        {0.90, 2.6366619772749E+01},  // time = 0.90s
        {0.95, 2.5948921136608E+01},  // time = 0.95s
        {1.00, 2.5190079591555E+01}   // time = 1.00s
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

    //-----------------------------------------------------
    // verify tip displacement history
    //-----------------------------------------------------

    // Using relative tolerance (10%) with minimum absolute threshold for small values
    // NOTE: The tolerance is high since we don't expect Kynema's constraint-based approach
    // to match BeamDyn's Dirichlet BCs based approach. Displacements are much lower than
    // the tolerance for majority of the time stamps.
    constexpr double relative_tolerance{0.10};  // 10% relative tolerance
    constexpr double min_abs_tolerance{1e-5};   // Minimum absolute tolerance for near-zero values
    auto expect_near_with_relative_tolerance = [&](const auto& actual, const auto& expected) {
        const auto tolerance = std::max(relative_tolerance * std::abs(expected), min_abs_tolerance);
        EXPECT_NEAR(actual, expected, tolerance);
    };

    // Verify all expected times are present and compare displacements
    for (const auto& [time, ref_disp] : beamdyn_displacements) {
        const auto it = tip_displacements.find(time);
        ASSERT_NE(it, tip_displacements.end()) << "Missing tip displacement data for time: " << time;
        expect_near_with_relative_tolerance(-it->second[2], ref_disp);
    }
}

}  // namespace kynema::tests
