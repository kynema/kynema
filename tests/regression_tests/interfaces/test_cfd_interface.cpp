#include <cmath>

#include <gtest/gtest.h>

#include "interfaces/cfd/interface.hpp"
#include "interfaces/cfd/interface_builder.hpp"

namespace kynema::tests {

using namespace kynema::interfaces::cfd;

TEST(CFDInterfaceTest, PrecessionTest) {
    // Create cfd interface
    constexpr auto mass_matrix =
        std::array{std::array{1., 0., 0., 0., 0., 0.}, std::array{0., 1., 0., 0., 0., 0.},
                   std::array{0., 0., 1., 0., 0., 0.}, std::array{0., 0., 0., 1., 0., 0.},
                   std::array{0., 0., 0., 0., 1., 0.}, std::array{0., 0., 0., 0., 0., .5}};
    const auto write_output{true};

    auto builder = InterfaceBuilder{}
                       .SetTimeStep(0.01)
                       .SetDampingFactor(1.)
                       .SetMaximumNonlinearIterations(5U)
                       .EnableFloatingPlatform(true)
                       .SetFloatingPlatformVelocity({0., 0., 0., 0.5, 0.5, 1.})
                       .SetFloatingPlatformMassMatrix(mass_matrix);

    if (write_output) {
        builder.SetOutputFile("CFDInterfaceTest.PrecessionTest");
    }

    auto interface = builder.Build();

    // Create reference to platform node in interface
    auto& platform_node = interface.turbine.floating_platform.node;

    // Run simulation for 500 steps
    for ([[maybe_unused]] auto i : std::views::iota(0, 500)) {
        EXPECT_EQ(interface.Step(), true);
    }

    // Check results at 500 steps
    EXPECT_NEAR(platform_node.displacement[0], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[1], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[2], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[3], -0.63053045128590757, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[4], 0.60556039120583116, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[5], -0.30157681970585326, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[6], -0.38049886257377241, 1.e-12);

    if (write_output) {
        // Read and verify data from NetCDF file
        const util::NetCdfFile file("CFDInterfaceTest.PrecessionTest/cfd_interface.nc", false);
        std::vector<double> displacements(1);

        // Check displacement at step 500, platform node
        const std::vector<size_t> start = {499, platform_node.id};
        const std::vector<size_t> count = {1, 1};

        file.ReadVariableAt("u_x", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], 0., 1.e-12);

        file.ReadVariableAt("u_y", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], 0., 1.e-12);

        file.ReadVariableAt("u_z", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], 0., 1.e-12);

        file.ReadVariableAt("u_w", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], -0.63053045128590757, 1.e-12);

        file.ReadVariableAt("u_i", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], 0.60556039120583116, 1.e-12);

        file.ReadVariableAt("u_j", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], -0.30157681970585326, 1.e-12);

        file.ReadVariableAt("u_k", start, count, displacements.data());
        EXPECT_NEAR(displacements[0], -0.38049886257377241, 1.e-12);
    }

    // Save the current state
    interface.SaveState();

    // Run simulation for an additional 100 steps
    for ([[maybe_unused]] auto iteration : std::views::iota(500, 600)) {
        EXPECT_EQ(interface.Step(), true);
    }

    // Check results at 600 steps
    EXPECT_NEAR(platform_node.displacement[0], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[1], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[2], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[3], -0.35839735533217154, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[4], 0.31963488733920042, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[5], -0.27587277763126217, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[6], -0.83263382410791731, 1.e-12);

    // Restore to saved state at 500 steps
    interface.RestoreState();

    // Check that results match 500 steps
    EXPECT_NEAR(platform_node.displacement[0], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[1], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[2], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[3], -0.63053045128590757, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[4], 0.60556039120583116, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[5], -0.30157681970585326, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[6], -0.38049886257377241, 1.e-12);

    // Run simulation from 500 to 600 steps
    for ([[maybe_unused]] auto interation : std::views::iota(500, 600)) {
        EXPECT_EQ(interface.Step(), true);
    }

    // Check that simulation gives same results at 600 steps
    EXPECT_NEAR(platform_node.displacement[0], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[1], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[2], 0., 1.e-12);
    EXPECT_NEAR(platform_node.displacement[3], -0.35839735533217154, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[4], 0.31963488733920042, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[5], -0.27587277763126217, 1.e-12);
    EXPECT_NEAR(platform_node.displacement[6], -0.83263382410791731, 1.e-12);
}

TEST(CFDInterfaceTest, FloatingPlatform) {
    // Solution parameters
    constexpr auto time_step = 0.01;
    constexpr auto t_end = 1.;
    constexpr auto rho_inf = 0.0;
    constexpr auto max_iter = 5;
    const auto n_steps{static_cast<size_t>(ceil(t_end / time_step)) + 1};
    const auto write_output{false};

    // Construct platform mass matrix
    constexpr auto platform_mass{1.419625E+7};                                     // kg
    constexpr auto gravity = std::array{0., 0., -9.8124};                          // m/s/s
    constexpr auto platform_moi = std::array{1.2898E+10, 1.2851E+10, 1.4189E+10};  // kg*m*m
    constexpr auto platform_mass_matrix = std::array{
        std::array{platform_mass, 0., 0., 0., 0., 0.},    // Row 1
        std::array{0., platform_mass, 0., 0., 0., 0.},    // Row 2
        std::array{0., 0., platform_mass, 0., 0., 0.},    // Row 3
        std::array{0., 0., 0., platform_moi[0], 0., 0.},  // Row 4
        std::array{0., 0., 0., 0., platform_moi[1], 0.},  // Row 5
        std::array{0., 0., 0., 0., 0., platform_moi[2]},  // Row 6
    };

    // Mooring line properties
    constexpr auto mooring_line_stiffness{48.9e3};       // N
    constexpr auto mooring_line_initial_length{55.432};  // m

    // Create cfd interface
    auto builder = InterfaceBuilder{}
                       .SetGravity(gravity)
                       .SetTimeStep(time_step)
                       .SetDampingFactor(rho_inf)
                       .SetMaximumNonlinearIterations(max_iter)
                       .EnableFloatingPlatform(true)
                       .SetFloatingPlatformPosition({0., 0., -7.53, 1., 0., 0., 0.})
                       .SetFloatingPlatformMassMatrix(platform_mass_matrix)
                       .SetNumberOfMooringLines(3)
                       .SetMooringLineStiffness(0, mooring_line_stiffness)
                       .SetMooringLineUndeformedLength(0, mooring_line_initial_length)
                       .SetMooringLineFairleadPosition(0, {-40.87, 0.0, -14.})
                       .SetMooringLineAnchorPosition(0, {-105.47, 0.0, -58.4})
                       .SetMooringLineStiffness(1, mooring_line_stiffness)
                       .SetMooringLineUndeformedLength(1, mooring_line_initial_length)
                       .SetMooringLineFairleadPosition(1, {20.43, -35.39, -14.})
                       .SetMooringLineAnchorPosition(1, {52.73, -91.34, -58.4})
                       .SetMooringLineStiffness(2, mooring_line_stiffness)
                       .SetMooringLineUndeformedLength(2, mooring_line_initial_length)
                       .SetMooringLineFairleadPosition(2, {20.43, 35.39, -14.})
                       .SetMooringLineAnchorPosition(2, {52.73, 91.34, -58.4});
    if (write_output) {
        builder.SetOutputFile("CFDInterfaceTest.FloatingPlatform");
    }

    auto interface = builder.Build();

    // Calculate buoyancy force as percentage of gravitational force plus spring forces times
    const auto initial_spring_force = 1907514.4912628897;
    const auto platform_gravity_force = -gravity[2] * platform_mass;
    const auto buoyancy_force = initial_spring_force + platform_gravity_force;

    // Iterate through time steps
    for (auto i : std::views::iota(0U, n_steps)) {
        // Calculate current time
        const auto t = static_cast<double>(i) * time_step;

        // Apply load in y direction
        interface.turbine.floating_platform.node.loads[1] =
            1e6 * sin(2. * std::numbers::pi / 20. * t);

        // Apply time varying buoyancy force
        interface.turbine.floating_platform.node.loads[2] =
            buoyancy_force + 0.5 * initial_spring_force * sin(2. * std::numbers::pi / 20. * t);

        // Apply time varying moments to platform node
        interface.turbine.floating_platform.node.loads[3] =
            5.0e5 * sin(2. * std::numbers::pi / 15. * t);  // rx
        interface.turbine.floating_platform.node.loads[4] =
            1.0e6 * sin(2. * std::numbers::pi / 30. * t);  // ry
        interface.turbine.floating_platform.node.loads[5] =
            2.0e7 * sin(2. * std::numbers::pi / 60. * t);  // rz

        // Step
        const auto converged = interface.Step();
        EXPECT_EQ(converged, true);

        // Check for expected displacements/rotations of platform node
        if (i == 100) {
            const auto& platform_u = interface.turbine.floating_platform.node.displacement;
            EXPECT_NEAR(platform_u[0], 3.4379940641493395e-06, 1e-12);
            EXPECT_NEAR(platform_u[1], 0.0036709089002912518, 1e-12);
            EXPECT_NEAR(platform_u[2], 0.0034998271929752248, 1e-12);
            EXPECT_NEAR(platform_u[3], 0.99999999992274302, 1e-12);
            EXPECT_NEAR(platform_u[4], 1.3448158662549561e-06, 1e-12);
            EXPECT_NEAR(platform_u[5], 1.2713890988278111e-06, 1e-12);
            EXPECT_NEAR(platform_u[6], 1.2291853012439073e-05, 1e-12);
        }
    }
}

TEST(CFDInterfaceTest, Restart) {
    // Solution parameters
    constexpr auto time_step = 0.01;
    constexpr auto rho_inf = 0.0;
    constexpr auto max_iter = 5;

    // Construct platform mass matrix
    constexpr auto platform_mass{1.419625E+7};                                     // kg
    constexpr auto gravity = std::array{0., 0., -9.8124};                          // m/s/s
    constexpr auto platform_moi = std::array{1.2898E+10, 1.2851E+10, 1.4189E+10};  // kg*m*m
    constexpr auto platform_mass_matrix = std::array{
        std::array{platform_mass, 0., 0., 0., 0., 0.},    // Row 1
        std::array{0., platform_mass, 0., 0., 0., 0.},    // Row 2
        std::array{0., 0., platform_mass, 0., 0., 0.},    // Row 3
        std::array{0., 0., 0., platform_moi[0], 0., 0.},  // Row 4
        std::array{0., 0., 0., 0., platform_moi[1], 0.},  // Row 5
        std::array{0., 0., 0., 0., 0., platform_moi[2]},  // Row 6
    };

    // Mooring line properties
    constexpr auto mooring_line_stiffness{48.9e3};       // N
    constexpr auto mooring_line_initial_length{55.432};  // m

    auto builder = InterfaceBuilder{}
                       .SetGravity(gravity)
                       .SetTimeStep(time_step)
                       .SetDampingFactor(rho_inf)
                       .SetMaximumNonlinearIterations(max_iter)
                       .EnableFloatingPlatform(true)
                       .SetFloatingPlatformPosition({0., 0., -7.53, 1., 0., 0., 0.})
                       .SetFloatingPlatformMassMatrix(platform_mass_matrix)
                       .SetNumberOfMooringLines(3)
                       .SetMooringLineStiffness(0, mooring_line_stiffness)
                       .SetMooringLineUndeformedLength(0, mooring_line_initial_length)
                       .SetMooringLineFairleadPosition(0, {-40.87, 0.0, -14.})
                       .SetMooringLineAnchorPosition(0, {-105.47, 0.0, -58.4})
                       .SetMooringLineStiffness(1, mooring_line_stiffness)
                       .SetMooringLineUndeformedLength(1, mooring_line_initial_length)
                       .SetMooringLineFairleadPosition(1, {20.43, -35.39, -14.})
                       .SetMooringLineAnchorPosition(1, {52.73, -91.34, -58.4})
                       .SetMooringLineStiffness(2, mooring_line_stiffness)
                       .SetMooringLineUndeformedLength(2, mooring_line_initial_length)
                       .SetMooringLineFairleadPosition(2, {20.43, 35.39, -14.})
                       .SetMooringLineAnchorPosition(2, {52.73, 91.34, -58.4});

    auto interface1 = builder.Build();

    // Take 10 initial steps
    for (auto i : std::views::iota(0U, 100U)) {
        const auto t = static_cast<double>(i) * time_step;
        interface1.turbine.floating_platform.node.loads[1] =
            1e6 * sin(2. * std::numbers::pi / 20. * t);
        auto converged = interface1.Step();
        EXPECT_TRUE(converged);
    }

    interface1.WriteRestart("test_restart.dat");

    // Take 10 more steps using original system
    for (auto i : std::views::iota(0U, 100U)) {
        const auto t = static_cast<double>(i) * time_step;
        interface1.turbine.floating_platform.node.loads[1] =
            1e6 * sin(2. * std::numbers::pi / 20. * t);
        auto converged = interface1.Step();
        EXPECT_TRUE(converged);
    }

    auto interface2 = builder.Build();
    interface2.ReadRestart("test_restart.dat");

    // Take 10 steps using restarted system
    for (auto i : std::views::iota(0U, 100U)) {
        const auto t = static_cast<double>(i) * time_step;
        interface2.turbine.floating_platform.node.loads[1] =
            1e6 * sin(2. * std::numbers::pi / 20. * t);
        auto converged = interface2.Step();
        EXPECT_TRUE(converged);
    }

    const auto& platform_u_1 = interface1.turbine.floating_platform.node.displacement;
    const auto& platform_u_2 = interface2.turbine.floating_platform.node.displacement;
    // Ensure platform location is the same for original and restarted system
    for (auto i : std::views::iota(0U, 7U)) {
        EXPECT_NEAR(platform_u_1[i], platform_u_2[i], 1.e-12);
    }

    std::filesystem::remove("test_restart.dat");
}

}  // namespace kynema::tests
