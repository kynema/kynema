#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

TEST(DynamicVerificationTest, Iea15MwBladeBending) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};
    const double time_step{0.01};
    builder.Solution()
        .EnableDynamicSolve()               // Dynamic analysis
        .SetTimeStep(time_step)             // Time step size
        .SetDampingFactor(0.)               // Max numerical damping (ρ_∞ = 0.)
        .SetMaximumNonlinearIterations(15)  // Max Newton-Raphson iterations
        .SetAbsoluteErrorTolerance(1e-7)    // Absolute error tolerance
        .SetRelativeErrorTolerance(1e-5);   // Relative error tolerance

    if (write_output) {
        builder.Outputs().SetOutputFilePath("DynamicVerificationTest.Iea15MwTurbineBladeBending");
    }

    //----------------------------------
    // beam element
    //----------------------------------
    const int num_nodes{15};                 // number of nodes = n
    const int element_order{num_nodes - 1};  // element order = n-1
    const int section_refinement{
        num_nodes - 1
    };  // each section of blade uses n-pt Gauss-Legendre quadrature for integration

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

    // Add reference axis twist (deg -> rad)
    const auto twist = wio_blade["outer_shape_bem"]["twist"];
    const auto twist_grid = twist["grid"].as<std::vector<double>>();
    const auto twist_values = twist["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, twist_grid.size())) {
        builder.Blade().AddRefAxisTwist(twist_grid[i], twist_values[i] * std::numbers::pi / 180.);
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

    //------------------------------------------------------
    // apply load and advance simulation for 10s
    //------------------------------------------------------
    auto& tip_node = interface.Blade().nodes.back();

    // Apply tip load in flapwise direction i.e. -y axis, enough to cause ~10% deflection of tip node
    tip_node.loads[1] = -2.e5;  // 200 kN

    /*
    std::cout << "Time, Tip node displacement in x direction, Tip node displacement in y direction, "
              << "Tip node displacement in z direction" << "\n";
    std::cout << 0. << ", " << std::setprecision(15) << tip_node.displacement[0] << ", "
              << tip_node.displacement[1] << ", " << tip_node.displacement[2] << "\n";
    */
    const auto num_steps = static_cast<size_t>(0.1 / time_step);  // 0.1 s at time step size = 0.01 s
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

    //-----------------------------------------------------
    // verify tip displacement history against SeaHOWL
    //-----------------------------------------------------

    // TODO
}

}  // namespace kynema::tests
