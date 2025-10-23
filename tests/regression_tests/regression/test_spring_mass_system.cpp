#include <numbers>

#include <gtest/gtest.h>

#include "model/model.hpp"
#include "step/step.hpp"
#include "step/step_parameters.hpp"

namespace {

/*
 * A simple spring-mass system with one mass attached to a fixed point by a spring.
 * The mass oscillates horizontally with simple harmonic motion. The system is initialized
 * with the mass displaced 2 units from equilibrium position.
 *
 * Schematic of the system with mass (M) and spring (/\/\/), fixed BC at left end node displacement
 * applied at right node:
 *
 * /|
 * /|o---/\/\/--- ( M )o---->
 * /|
 */
template <typename DeviceType>
inline auto SetUpSpringMassSystem() {
    auto model = kynema::Model();

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

    // Add fixed BC to the first node
    model.AddFixedBC(fixed_node_id);

    // Add spring element
    const auto k = 10.;  // stiffness
    const auto l0 = 0.;  // undeformed length
    model.AddSpringElement(fixed_node_id, mass_node_id, k, l0);

    // The spring-mass system should move periodically between -2 and 2 (position of the second node)
    // with simple harmonic motion where the time period is 2 * pi * sqrt(m / k)
    const double T = 2. * std::numbers::pi * sqrt(m / k);
    constexpr auto num_steps = 1000;

    // Set up step parameters
    constexpr bool is_dynamic_solve(true);
    constexpr size_t max_iter(6);
    constexpr double rho_inf(0.);                                // No damping
    const double step_size(T / static_cast<double>(num_steps));  // Calculate step size
    auto parameters = kynema::StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints, solver] = model.CreateSystemWithSolver<DeviceType>();

    // Create host mirror for checking solution
    auto q = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.q);

    // Run simulation for T seconds
    for (auto time_step : std::views::iota(1U, num_steps + 1U)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);
        // Simulation at time T / 2
        if (time_step == num_steps / 2) {
            Kokkos::deep_copy(q, state.q);
            EXPECT_NEAR(q(0, 0), 0., 1.e-14);  // First node is fixed
            EXPECT_NEAR(
                q(1, 0), -3.9999199193098396, 1.e-12
            );  // Second node should have displacement close to -4.0 after T/2
        }
        // Simulation at time T
        if (time_step == num_steps) {
            Kokkos::deep_copy(q, state.q);
            EXPECT_NEAR(q(0, 0), 0., 1.e-14);  // First node is fixed
            EXPECT_NEAR(
                q(1, 0), -8.1226588438437419e-05, 1.e-12
            );  // Second node should have displacement close to 0. after T
        }
    }
}

/*
 * A chain of identical masses held together by identical springs and anchored at both ends.
 * When the masses are equidistant from eachother, the system shoud stay in perfect static
 * equilibrium. By adjusting the variable number_of_masses, the size of the system can be scaled up
 * or down.
 *
 * ASCII Art of the system with three masses (M) and anchor points (A):
 *
 * A ---/\/\/--- M ---/\/\/--- M ---/\/\/--- M ---/\/\/--- A
 */
inline auto SetUpSpringMassChainSystem() {
    auto model = kynema::Model();

    // Add nodes for each mass and an anchor point on each side
    constexpr auto number_of_masses = 10U;
    constexpr auto displacement = 0.5;
    auto position = 0.;
    model.AddNode().SetPosition(position, 0., 0., 1., 0., 0., 0.);
    for ([[maybe_unused]] auto mass_number : std::views::iota(0U, number_of_masses)) {
        position += displacement;
        model.AddNode().SetPosition(position, 0., 0., 1., 0., 0., 0.);
    }
    position += displacement;
    model.AddNode().SetPosition(position, 0., 0., 1., 0., 0., 0.);

    // Mass matrix (Identical for all masses)
    constexpr auto m = 1.;
    constexpr auto j = 1.;
    constexpr auto mass_matrix = std::array{
        std::array{m, 0., 0., 0., 0., 0.},  // mass in x-direction
        std::array{0., m, 0., 0., 0., 0.},  // mass in y-direction
        std::array{0., 0., m, 0., 0., 0.},  // mass in z-direction
        std::array{0., 0., 0., j, 0., 0.},  // inertia around x-axis
        std::array{0., 0., 0., 0., j, 0.},  // inertia around y-axis
        std::array{0., 0., 0., 0., 0., j},  // inertia around z-axis
    };

    for (auto mass_number : std::views::iota(0U, number_of_masses)) {
        model.AddMassElement(mass_number + 1, mass_matrix);
    }

    // Create springs
    const auto k = 10.;
    for (auto mass_number : std::views::iota(0U, number_of_masses + 1U)) {
        model.AddSpringElement(mass_number, mass_number + 1, k, 0.);
    }

    // Add fixed BC to the anchor nodes
    model.AddFixedBC(0);
    model.AddFixedBC(number_of_masses + 1);

    // Create system and solver
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Set up step parameters
    constexpr auto num_steps = 1000;
    constexpr bool is_dynamic_solve(true);
    constexpr size_t max_iter(6);
    constexpr double rho_inf(0.);
    const double final_time = 2. * std::numbers::pi * sqrt(m / k);
    const double step_size(final_time / static_cast<double>(num_steps));
    auto parameters = kynema::StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    auto q = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.q);

    // Run simulation for T seconds
    for ([[maybe_unused]] auto time_step : std::views::iota(1U, num_steps + 1U)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);
        Kokkos::deep_copy(q, state.q);
        for (auto node : std::views::iota(0U, q.extent(0))) {
            EXPECT_NEAR(q(node, 0), 0., 1.e-14);
        }
    }
}

}  // namespace

namespace kynema::tests {

TEST(SpringMassSystemTest, FinalDisplacement) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    SetUpSpringMassSystem<DeviceType>();
}

TEST(SpringMassSystemTest, FinalDisplacement_Host) {
    using DeviceType = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
    SetUpSpringMassSystem<DeviceType>();
}

TEST(SpringMassChainSystemTest, FinalDisplacement) {
    SetUpSpringMassChainSystem();
}

}  // namespace kynema::tests
