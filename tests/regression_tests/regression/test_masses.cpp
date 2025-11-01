#include <gtest/gtest.h>

#include "model/model.hpp"
#include "state/state.hpp"
#include "step/step.hpp"
#include "test_utilities.hpp"

namespace {
inline auto SetUpMasses() {
    // Create model object
    auto model = kynema::Model();

    // Set gravity
    model.SetGravity(0., 0., 9.81);

    // Add node
    const auto node_id = model.AddNode().SetPosition(0, 0, 0, 1, 0, 0, 0).Build();

    // Add mass element
    model.AddMassElement(
        node_id, {{
                     {1., 0., 0., 0., 0., 0.},
                     {0., 1., 0., 0., 0., 0.},
                     {0., 0., 1., 0., 0., 0.},
                     {0., 0., 0., 1., 0., 0.},
                     {0., 0., 0., 0., 1., 0.},
                     {0., 0., 0., 0., 0., 1.},
                 }}
    );

    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    // Initialize masses
    auto masses = model.CreateMasses<DeviceType>();

    // Create state
    auto state = model.CreateState<DeviceType>();

    auto parameters = kynema::StepParameters(false, 0, 0., 0.);
    kynema::step::UpdateSystemVariablesMasses(parameters, masses, state);

    return masses;
}

}  // namespace

namespace kynema::tests {

TEST(MassesTest, NodeInitialPosition) {
    const auto masses = SetUpMasses();
    expect_kokkos_view_2D_equal(
        masses.node_x0,
        {
            {0., 0., 0., 1., 0., 0., 0.},
        }
    );
}

TEST(MassesTest, MassMatrixInMaterialFrame) {
    const auto masses = SetUpMasses();
    expect_kokkos_view_2D_equal(
        Kokkos::subview(masses.qp_Mstar, 0, Kokkos::ALL, Kokkos::ALL),
        {
            {1., 0., 0., 0., 0., 0.},
            {0., 1., 0., 0., 0., 0.},
            {0., 0., 1., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0.},
            {0., 0., 0., 0., 1., 0.},
            {0., 0., 0., 0., 0., 1.},
        }
    );
}

TEST(MassesTest, GravityVector) {
    const auto masses = SetUpMasses();
    expect_kokkos_view_1D_equal(masses.gravity, {0., 0., 9.81});
}

TEST(MassesTest, ExternalForce) {
    // Create model object
    auto model = Model();

    // Set gravity
    model.SetGravity(0., 0., 0.);

    // Add node
    const auto node_id = model.AddNode().SetPosition(0, 0, 0, 1, 0, 0, 0).Build();

    // Add mass element
    const auto m = 1.;
    const auto j = 1.;
    model.AddMassElement(
        node_id, {{
                     {m, 0., 0., 0., 0., 0.},
                     {0., m, 0., 0., 0., 0.},
                     {0., 0., m, 0., 0., 0.},
                     {0., 0., 0., j, 0., 0.},
                     {0., 0., 0., 0., j, 0.},
                     {0., 0., 0., 0., 0., j},
                 }}
    );

    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    // Initialize masses
    auto masses = model.CreateMasses<DeviceType>();

    // Create state
    auto state = model.CreateState<DeviceType>();

    // Create solution parameters
    const auto time_step = 0.001;
    auto parameters = StepParameters(true, 5, time_step, 0.0);

    auto constraints = model.CreateConstraints<DeviceType>();
    auto elements = model.CreateElements<
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>>(
    );

    auto solver = CreateSolver<>(state, elements, constraints);

    const auto force_x = 5.;
    const auto torque_y = 3.;
    auto host_f = Kokkos::create_mirror_view(state.f);
    host_f(node_id, 0) = force_x;
    host_f(node_id, 4) = torque_y;
    Kokkos::deep_copy(state.f, host_f);

    // Run simulation for 1000 steps
    const auto n_steps = 1000;
    for ([[maybe_unused]] auto iteration : std::views::iota(0, n_steps)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);
    }

    // Get velocity
    auto v_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.v);

    EXPECT_NEAR(v_host(0, 0), 4.9975, 1.e-12);
    EXPECT_NEAR(v_host(0, 1), 0., 1.e-12);
    EXPECT_NEAR(v_host(0, 2), 0., 1.e-12);
    EXPECT_NEAR(v_host(0, 3), 0., 1.e-12);
    EXPECT_NEAR(v_host(0, 4), 2.9985, 1.e-12);
    EXPECT_NEAR(v_host(0, 5), 0., 1.e-12);
}

}  // namespace kynema::tests
