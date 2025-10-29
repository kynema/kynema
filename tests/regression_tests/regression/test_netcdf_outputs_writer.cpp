#include <filesystem>
#include <numbers>

#include <gtest/gtest.h>

#include "model/model.hpp"
#include "step/step.hpp"
#include "utilities/netcdf/node_state_writer.hpp"

namespace kynema::tests {

TEST(NetCDFOutputsWriterTest, SpringMassSystemOutputs) {
    auto model = Model();

    // Add two nodes for the spring element
    const auto fixed_node_id =
        model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();  // First node at origin
    const auto mass_node_id =
        model.AddNode().SetPosition(2., 0., 0., 1., 0., 0., 0.).Build();  // Second node at (2,0,0)

    // Add mass element
    constexpr auto m = 1.;
    constexpr auto j = 1.;
    model.AddMassElement(
        mass_node_id, {{
                          {m, 0., 0., 0., 0., 0.},  // mass in x-direction
                          {0., m, 0., 0., 0., 0.},  // mass in y-direction
                          {0., 0., m, 0., 0., 0.},  // mass in z-direction
                          {0., 0., 0., j, 0., 0.},  // inertia around x-axis
                          {0., 0., 0., 0., j, 0.},  // inertia around y-axis
                          {0., 0., 0., 0., 0., j},  // inertia around z-axis
                      }}
    );

    // Add fixed BC and spring element
    model.AddFixedBC(fixed_node_id);
    const auto k = 10.;  // stiffness
    const auto l0 = 0.;  // undeformed length
    model.AddSpringElement(fixed_node_id, mass_node_id, k, l0);

    // Set up simulation parameters
    const double T = 2. * std::numbers::pi * sqrt(m / k);
    constexpr auto num_steps = 1000;
    constexpr bool is_dynamic_solve(true);
    constexpr size_t max_iter(6);
    constexpr double rho_inf(0.);
    const double step_size(T / static_cast<double>(num_steps));
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver and system
    auto [state, elements, constraints, solver] = model.CreateSystemWithSolver();
    auto q = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.q);

    // Set up NodeStateWriter with 2 nodes
    const std::string output_file = "spring_mass_system.nc";
    std::filesystem::remove(output_file);  // Remove any existing file
    util::NodeStateWriter writer(
        output_file,
        true,   // create new file
        2,      // 2 nodes
        {"u"},  // write only displacements
        false,  // do not write deformation
        0       // no buffering
    );

    // Run simulation and write outputs
    for (auto time_step : std::views::iota(0U, num_steps + 1U)) {
        Kokkos::deep_copy(q, state.q);

        // Write only displacement data
        writer.WriteStateDataAtTimestep(
            time_step,
            "u",                 // displacement
            {q(0, 0), q(1, 0)},  // x component
            {q(0, 1), q(1, 1)},  // Empty
            {q(0, 2), q(1, 2)},  // Empty
            {q(0, 3), q(1, 3)},  // Empty
            {q(0, 4), q(1, 4)},  // Empty
            {q(0, 5), q(1, 5)},  // Empty
            {q(0, 6), q(1, 6)}   // Empty
        );

        // Step the simulation forward
        if (time_step < num_steps) {
            auto converged = Step(parameters, solver, elements, state, constraints);
            EXPECT_TRUE(converged);
        }
    }

    EXPECT_TRUE(std::filesystem::exists(output_file));

    // Verify netcdf output data
    const util::NetCdfFile file(output_file, false);
    std::vector<double> x_displacements(2);

    // Check displacement at T/2
    const std::vector<size_t> start_t_half = {num_steps / 2, 0};
    const std::vector<size_t> count = {1, 2};
    file.ReadVariableAt("u_x", start_t_half, count, x_displacements.data());
    EXPECT_NEAR(x_displacements[0], 0., 1.e-14);                   // First node is fixed
    EXPECT_NEAR(x_displacements[1], -3.9999199193098396, 1.e-12);  // Second node at T/2

    // Check displacement at T
    const std::vector<size_t> start_t = {num_steps, 0};
    file.ReadVariableAt("u_x", start_t, count, x_displacements.data());
    EXPECT_NEAR(x_displacements[0], 0., 1.e-14);  // First node is fixed
    EXPECT_NEAR(x_displacements[1], -8.1226588438437419e-05, 1.e-12);

    std::filesystem::remove(output_file);
}

}  // namespace kynema::tests
