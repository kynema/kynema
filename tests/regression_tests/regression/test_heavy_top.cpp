#include <gtest/gtest.h>

#include "model/model.hpp"
#include "step/step.hpp"
#include "test_utilities.hpp"

namespace kynema::tests {

inline auto SetUpHeavyTopTest() {
    auto model = Model();
    model.SetGravity(0., 0., -9.81);

    // Heavy top model parameters
    constexpr auto mass = 15.;                                           // mass
    constexpr auto inertia = std::array{0.234375, 0.46875, 0.234375};    // inertia matrix
    const auto x = Eigen::Matrix<double, 3, 1>(0., 1., 0.);              // initial position
    const auto omega = Eigen::Matrix<double, 3, 1>(0., 150., -4.61538);  // initial angular velocity
    const auto x_dot = omega.cross(x);                                   // initial velocity
    const auto omega_dot =
        std::array{661.3461692307691919, 0., 0.};  // initial anguluar acceleration
    const auto x_ddot =
        std::array{0., -21.3017325444000001, -30.9608307692308244};  // initial acceleration

    EXPECT_NEAR(x_dot[0], 4.61538, 1.e-15);
    EXPECT_NEAR(x_dot[1], 0., 1.e-15);
    EXPECT_NEAR(x_dot[2], 0., 1.e-15);

    // Add node with initial position and velocity
    auto mass_node_id =
        model.AddNode()
            .SetPosition(x(0), x(1), x(2), 1., 0., 0., 0.)
            .SetVelocity(x_dot(0), x_dot(1), x_dot(2), omega(0), omega(1), omega(2))
            .SetAcceleration(
                x_ddot[0], x_ddot[1], x_ddot[2], omega_dot[0], omega_dot[1], omega_dot[2]
            )
            .Build();

    // Add masses element with m and J as mass matrix
    model.AddMassElement(
        mass_node_id, {{
                          {mass, 0., 0., 0., 0., 0.},        // mass in x-direction
                          {0., mass, 0., 0., 0., 0.},        // mass in y-direction
                          {0., 0., mass, 0., 0., 0.},        // mass in z-direction
                          {0., 0., 0., inertia[0], 0., 0.},  // inertia xx
                          {0., 0., 0., 0., inertia[1], 0.},  // inertia yy
                          {0., 0., 0., 0., 0., inertia[2]},  // inertia zz
                      }}
    );

    // Add ground node at origin
    auto ground_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();

    // Add constraints (6 DOF base node -> 3 DOF target node)
    model.AddRigidJoint6DOFsTo3DOFs(std::array{mass_node_id, ground_node_id});
    model.AddPrescribedBC3DOFs(ground_node_id);

    // Set up step parameters
    constexpr bool is_dynamic_solve(true);
    constexpr size_t max_iter(10);
    constexpr double step_size(0.002);
    constexpr double rho_inf(0.9);
    constexpr double a_tol(1e-5);
    constexpr double r_tol(1e-3);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf, a_tol, r_tol);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints, solver] = model.CreateSystemWithSolver<>();

    // Run simulation for 0.8 seconds
    for ([[maybe_unused]] auto i : std::views::iota(0, 400)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);
    }

    const auto q = kokkos_view_2D_to_vector(state.q);

    EXPECT_NEAR(q[0][0], -0.42217802273894345, 1e-10);
    EXPECT_NEAR(q[0][1], -0.09458263530050703, 1e-10);
    EXPECT_NEAR(q[0][2], -0.04455460488952848, 1e-10);
    EXPECT_NEAR(q[0][3], -0.17919607435565366, 1e-10);
    EXPECT_NEAR(q[0][4], 0.21677896640311572, 1e-10);
    EXPECT_NEAR(q[0][5], -0.9594776960853596, 1e-10);
    EXPECT_NEAR(q[0][6], -0.017268392381761217, 1e-10);
}

TEST(HeavyTopTest, FinalState) {
    SetUpHeavyTopTest();
}

}  // namespace kynema::tests
