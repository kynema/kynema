#include <gtest/gtest.h>

#include "model/model.hpp"
#include "step/step.hpp"

namespace {
inline auto SetUpPrecessionTest() {
    auto model = kynema::Model();

    // Add node with initial position and velocity
    auto node_id = model.AddNode().SetVelocity(0., 0., 0., 0.5, 0.5, 1.0).Build();

    // Add masses element
    constexpr auto m = 1.0;
    model.AddMassElement(
        node_id, {{
                     {m, 0., 0., 0., 0., 0.},    //
                     {0., m, 0., 0., 0., 0.},    //
                     {0., 0., m, 0., 0., 0.},    //
                     {0., 0., 0., 1., 0., 0.},   //
                     {0., 0., 0., 0., 1., 0.},   //
                     {0., 0., 0., 0., 0., 0.5},  //
                 }}
    );

    // Set up step parameters
    constexpr bool is_dynamic_solve(true);
    constexpr size_t max_iter(6);
    constexpr double step_size(0.01);
    constexpr double rho_inf(1.0);
    auto parameters = kynema::StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints, solver] = model.CreateSystemWithSolver<>();

    // Run simulation for 500 steps
    for ([[maybe_unused]] auto i : std::views::iota(0, 500)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);
    }

    auto q_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.q);
    EXPECT_NEAR(q_host(0, 0), 0., 1.e-12);
    EXPECT_NEAR(q_host(0, 1), 0., 1.e-12);
    EXPECT_NEAR(q_host(0, 2), 0., 1.e-12);
    EXPECT_NEAR(q_host(0, 3), -0.63053045128590757, 1.e-12);
    EXPECT_NEAR(q_host(0, 4), 0.60556039120583116, 1.e-12);
    EXPECT_NEAR(q_host(0, 5), -0.30157681970585326, 1.e-12);
    EXPECT_NEAR(q_host(0, 6), -0.38049886257377241, 1.e-12);
}

}  // namespace

namespace kynema::tests {

TEST(PrecessionTest, FinalRotation) {
    SetUpPrecessionTest();
}

}  // namespace kynema::tests
