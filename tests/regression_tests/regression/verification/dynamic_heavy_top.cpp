#include <numbers>

#include <gtest/gtest.h>

#include "model/model.hpp"
#include "step/step.hpp"

namespace kynema::tests {

TEST(VerificationTest, Dynamic_HeavyTopSpinningUnderGravity) {
    //----------------------------------
    // solution parameters
    //----------------------------------
    constexpr bool is_dynamic_solve{true};  // Dynamic analysis
    constexpr size_t max_iter{15};          // Max Newton-Raphson iterations
    constexpr double time_step{2e-4};       // Time step size (s)
    constexpr double rho_inf{0.9};          // Light numerical damping (ρ_∞ = 0.9)
    constexpr double a_tol{1e-7};           // Absolute error tolerance
    constexpr double r_tol{1e-6};           // Relative error tolerance
    auto parameters = StepParameters(is_dynamic_solve, max_iter, time_step, rho_inf, a_tol, r_tol);

    //----------------------------------
    // heavy top model parameters
    //----------------------------------
    // Heavy top is a rotating rigid body fixed to ground by a spherical joint
    constexpr auto mass = 15.;
    constexpr auto inertia = std::array{0.234375, 0.46875, 0.234375};

    // Initial conditions
    const auto x = Eigen::Matrix<double, 3, 1>(0., 1., 0.);              // Initial position
    const auto omega = Eigen::Matrix<double, 3, 1>(0., 150., -4.61538);  // Initial angular velocity
    const auto x_dot = omega.cross(x);                                   // Initial velocity

    // Initial accelerations (converted to inertial frame)
    const auto omega_dot_inertial_frame = std::array{-692.307, 0., 0.};
    const auto x_ddot_inertial_frame = std::array{0., 0., -0.654};

    //----------------------------------
    // build model
    //----------------------------------
    auto model = Model();
    model.SetGravity(0., 0., -9.81);

    // Add mass node with initial position, velocity, and acceleration
    auto mass_node_id =
        model.AddNode()
            .SetPosition(x[0], x[1], x[2], 1., 0., 0., 0.)
            .SetVelocity(x_dot[0], x_dot[1], x_dot[2], omega[0], omega[1], omega[2])
            .SetAcceleration(
                x_ddot_inertial_frame[0], x_ddot_inertial_frame[1], x_ddot_inertial_frame[2],
                omega_dot_inertial_frame[0], omega_dot_inertial_frame[1], omega_dot_inertial_frame[2]
            )
            .Build();

    // Add mass element with 6x6 mass matrix
    model.AddMassElement(
        mass_node_id, {{
                          {mass, 0., 0., 0., 0., 0.},        //
                          {0., mass, 0., 0., 0., 0.},        //
                          {0., 0., mass, 0., 0., 0.},        //
                          {0., 0., 0., inertia[0], 0., 0.},  //
                          {0., 0., 0., 0., inertia[1], 0.},  //
                          {0., 0., 0., 0., 0., inertia[2]},  //
                      }}
    );

    // Add ground node at origin (fixed)
    auto ground_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();

    // Add constraints: rigid joint (spherical joint) connecting mass to ground
    // 6 DOF base node (mass) -> 3 DOF target node (ground)
    model.AddRigidJoint6DOFsTo3DOFs(std::array{mass_node_id, ground_node_id});
    model.AddPrescribedBC3DOFs(ground_node_id);

    // Create system with solver
    auto [state, elements, constraints, solver] = model.CreateSystemWithSolver<>();

    //------------------------------------------------------
    // run simulation for 2 seconds
    //------------------------------------------------------
    constexpr double simulation_time{2.};  // s
    const auto num_steps = static_cast<size_t>(std::round(simulation_time / time_step));
    for ([[maybe_unused]] auto step : std::views::iota(1U, num_steps + 1U)) {
        auto converged = Step(parameters, solver, elements, state, constraints);

        // Verify we reach convergence
        ASSERT_EQ(converged, true);
    }

    //-----------------------------------------------------
    // verify results
    //-----------------------------------------------------

    // Benchmark solution calculated with Brüls, Cardona, and Arnold (2012) implementation
    // in material coordinate system with dt = 2e-6 s
    // Reference results from BeamDyn simulation stored in: heavy_top_kynema_material_0pt000002s.csv
    const std::array<double, 7> benchmark_solution_at_2_seconds{
        -0.150416985652933,  // u0
        -1.724164450779749,  // u1
        -0.673023460701136,  // u2
        0.122947275600560,   // q0
        -0.491079359219981,  // q1
        0.350430794926147,   // q2
        -0.787986857972906   // q3
    };

    // Using absolute tolerance based on observed numerical errors
    // Maximum observed error is ~2.05e-4, rounding up to 5.e-4
    constexpr double tolerance{5e-4};

    auto qHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.q);
    EXPECT_NEAR(qHost(mass_node_id, 0), benchmark_solution_at_2_seconds[0], tolerance);
    EXPECT_NEAR(qHost(mass_node_id, 1), benchmark_solution_at_2_seconds[1], tolerance);
    EXPECT_NEAR(qHost(mass_node_id, 2), benchmark_solution_at_2_seconds[2], tolerance);
    EXPECT_NEAR(qHost(mass_node_id, 3), benchmark_solution_at_2_seconds[3], tolerance);
    EXPECT_NEAR(qHost(mass_node_id, 4), benchmark_solution_at_2_seconds[4], tolerance);
    EXPECT_NEAR(qHost(mass_node_id, 5), benchmark_solution_at_2_seconds[5], tolerance);
    EXPECT_NEAR(qHost(mass_node_id, 6), benchmark_solution_at_2_seconds[6], tolerance);
}

}  // namespace kynema::tests
