#include <fstream>
#include <initializer_list>
#include <iostream>
#include <numbers>

#include <gtest/gtest.h>

#include "model/model.hpp"
#include "step/step.hpp"
#include "test_utilities.hpp"

namespace kynema::tests {

template <typename T>
void WriteMatrixToFile(const std::vector<std::vector<T>>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << "\n";
        return;
    }
    for (const auto& innerVector : data) {
        for (const auto& element : innerVector) {
            file << element << ",";
        }
        file << "\n";
    }
    file.close();
}

// Mass matrix for uniform composite beam section
constexpr auto mass_matrix = std::array{
    std::array{8.538e-2, 0., 0., 0., 0., 0.},    //
    std::array{0., 8.538e-2, 0., 0., 0., 0.},    //
    std::array{0., 0., 8.538e-2, 0., 0., 0.},    //
    std::array{0., 0., 0., 1.4433e-2, 0., 0.},   //
    std::array{0., 0., 0., 0., 0.40972e-2, 0.},  //
    std::array{0., 0., 0., 0., 0., 1.0336e-2},
};

// create a unity mass matrix
constexpr auto mass_matrix_unity = std::array{
    std::array{1., 0., 0., 0., 0., 0.},  //
    std::array{0., 1., 0., 0., 0., 0.},  //
    std::array{0., 0., 1., 0., 0., 0.},  //
    std::array{0., 0., 0., 1., 0., 0.},  //
    std::array{0., 0., 0., 0., 1., 0.},  //
    std::array{0., 0., 0., 0., 0., 1.},
};

// Stiffness matrix for uniform composite beam section
constexpr auto stiffness_matrix = std::array{
    std::array{1368.17e3, 0., 0., 0., 0., 0.},
    std::array{0., 88.56e3, 0., 0., 0., 0.},
    std::array{0., 0., 38.78e3, 0., 0., 0.},
    std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
    std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
    std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
};

// create a unit stiffness matrix
constexpr auto stiffness_matrix_unity = std::array{
    std::array{1., 0., 0., 0., 0., 0.},  //
    std::array{0., 1., 0., 0., 0., 0.},  //
    std::array{0., 0., 1., 0., 0., 0.},  //
    std::array{0., 0., 0., 1., 0., 0.},  //
    std::array{0., 0., 0., 0., 1., 0.},  //
    std::array{0., 0., 0., 0., 0., 1.},
};

// Node locations (GLL quadrature) - 6 nodes
const auto node_s = std::vector{
    0., 0.11747233803526763, 0.35738424175967748, 0.64261575824032247, 0.88252766196473242, 1.
};

// Element quadrature
const auto quadrature = std::vector<std::array<double, 2>>{
    {-0.9491079123427585, 0.1294849661688697},  {-0.7415311855993943, 0.27970539148927664},
    {-0.40584515137739696, 0.3818300505051189}, {6.123233995736766e-17, 0.4179591836734694},
    {0.4058451513773971, 0.3818300505051189},   {0.7415311855993945, 0.27970539148927664},
    {0.9491079123427585, 0.1294849661688697},
};

const auto sections = std::vector{
    BeamSection(0., mass_matrix, stiffness_matrix),
    BeamSection(1., mass_matrix, stiffness_matrix),
};

const auto sections_unity = std::vector{
    BeamSection(0., mass_matrix_unity, stiffness_matrix_unity),
    BeamSection(1., mass_matrix_unity, stiffness_matrix_unity),
};

TEST(RotatingBeamTest, StepConvergence) {
    auto model = Model();

    // Set gravity
    model.SetGravity(0., 0., 0.);

    const auto x0_root = Eigen::Matrix<double, 3, 1>(2., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    const double omega = 0.1;
    std::vector<size_t> node_ids;
    std::ranges::transform(node_s, std::back_inserter(node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode().SetElemLocation(s).SetPosition(x, 0., 0., 1., 0., 0., 0.).Build();
    });

    // Add beam element and set its position and velocity
    model.AddBeamElement(node_ids, sections, quadrature);
    model.TranslateBeam(0, std::array{x0_root(0), 0., 0.});
    model.SetBeamVelocityAboutPoint(
        0, std::array{0., 0., 0., 0., 0., omega}, std::array{0., 0., 0.}
    );

    // Add prescribed boundary condition on first node of beam
    model.AddPrescribedBC(node_ids[0]);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Create solver parameters
    const bool is_dynamic_solve(true);
    const int max_iter(5);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Perform 10 time steps and check for convergence within max_iter iterations
    for (auto i : std::views::iota(0, 10)) {
        // Set constraint displacement
        const auto q = Eigen::Quaternion<double>(Eigen::AngleAxis<double>(
            omega * step_size * (i + 1), Eigen::Matrix<double, 3, 1>::Unit(2)
        ));
        const auto x_root = q._transformVector(x0_root);
        const auto u_root = x_root - x0_root;
        const auto displacement =
            std::array{u_root(0), u_root(1), u_root(2), q.w(), q.x(), q.y(), q.z()};
        constraints.UpdateDisplacement(0, displacement);
        const auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_EQ(converged, true);
    }

    expect_kokkos_view_2D_equal(
        state.q,
        {
            {-0.000099999166669473282, 0.019999666668333329, 0., 0.99998750002604219, 0., 0.,
             0.0049999791666927116},
            {-0.00015873054129883523, 0.031746705307316693, 5.3452180213890101E-13,
             0.99998750002475689, 6.9132808328819882E-12, -1.2067607167192481E-13,
             0.0049999794236664759},
            {-0.00027867990247669591, 0.055737498222589707, 2.7206029385761709E-12,
             0.99998750002344106, 1.6236092075463668E-11, -1.0204079036815215E-13,
             0.004999979686833151},
            {-0.00042128940140729969, 0.084260177944726039, 2.2262966672810607E-12,
             0.99998750002292091, 2.3329848595908065E-11, 8.6195395416293551E-13,
             0.0049999797908453587},
            {-0.0005412402267101781, 0.10825097167135142, -5.7646870281175913E-12,
             0.99998750002269443, 2.5986319594619517E-11, 1.8229927487844355E-12,
             0.0049999798361397927},
            {-0.0005999751384186499, 0.11999801130335846, -1.0588553074456064E-11,
             0.99998750002267489, 2.6230990329572291E-11, 1.8800830187025051E-12,
             0.0049999798400182259},
        }
    );
}

inline void CreateTwoBeamSolverWithSameBeamsAndStep() {
    // Create model for managing nodes and constraints
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    constexpr auto num_blades = 2;
    constexpr auto velocity = std::array{0., 0., 0., 0., 0., 1.};
    constexpr auto origin = std::array{0., 0., 0.};
    constexpr auto hub_radius = 2.;
    for ([[maybe_unused]] auto blade_number : std::views::iota(0, num_blades)) {
        auto beam_node_ids = std::vector<size_t>(node_s.size());
        std::ranges::transform(node_s, std::begin(beam_node_ids), [&](auto s) {
            return model.AddNode()
                .SetElemLocation(s)
                .SetPosition(10. * s, 0., 0., 1., 0., 0., 0.)
                .Build();
        });
        auto blade_elem_id = model.AddBeamElement(beam_node_ids, sections, quadrature);
        auto rotation_quaternion = std::array{1., 0., 0., 0.};
        model.TranslateBeam(blade_elem_id, std::array{hub_radius, 0., 0.});
        model.RotateBeamAboutPoint(blade_elem_id, rotation_quaternion, origin);
        model.SetBeamVelocityAboutPoint(blade_elem_id, velocity, origin);
    }

    // Add a prescribed BC for each root node
    for (const auto& beam_element : model.GetBeamElements()) {
        model.AddPrescribedBC(beam_element.node_ids.front());
    }

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(1);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Calculate hub rotation for this time step
    const auto q_hub = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(step_size, Eigen::Matrix<double, 3, 1>(&velocity[3]))
    );

    // Define hub translation/rotation displacement
    const auto u_hub = std::array{0., 0., 0., q_hub.w(), q_hub.x(), q_hub.y(), q_hub.z()};

    // Update constraint displacements
    for (auto j : std::views::iota(0U, constraints.num_constraints)) {
        constraints.UpdateDisplacement(j, u_hub);
    }

    // Take step, don't check for convergence, the following tests check that
    // all the elements were assembled properly
    Step(parameters, solver, elements, state, constraints);

    auto n = solver.num_system_dofs / 2;
    auto m = constraints.num_dofs / 2;

    // Check that R vector is the same for both beams
    auto b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), solver.b);
    for (auto i : std::views::iota(0U, n)) {
        EXPECT_NEAR(b(i, 0), b(n + i, 0), 1.e-10);
    }

    // Check that Phi vector is the same for both beams
    auto Phi = kokkos_view_2D_to_vector(constraints.residual_terms);
    for (auto i : std::views::iota(0U, m)) {
        EXPECT_NEAR(Phi[0][i], Phi[1][i], 1.e-10);
    }
}

TEST(RotatingBeamTest, TwoBeam) {
    CreateTwoBeamSolverWithSameBeamsAndStep();
}

TEST(RotatingBeamTest, ThreeBladeRotor) {
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 9.81);

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 1 rad/s angular velocity around the z axis
    constexpr auto num_blades = 3;
    constexpr auto velocity = std::array{0., 0., 0., 0., 0., 1.};
    constexpr auto origin = std::array{0., 0., 0.};
    constexpr auto hub_radius = 2.;
    for (auto blade_number : std::views::iota(0, num_blades)) {
        auto beam_node_ids = std::vector<size_t>(node_s.size());
        std::ranges::transform(node_s, std::begin(beam_node_ids), [&](auto s) {
            return model.AddNode()
                .SetElemLocation(s)
                .SetPosition(10. * s, 0., 0., 1., 0., 0., 0.)
                .Build();
        });
        auto blade_elem_id = model.AddBeamElement(beam_node_ids, sections, quadrature);
        const auto blade_angle = 2. * std::numbers::pi * blade_number / num_blades;
        const auto rotation_quaternion = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(blade_angle, Eigen::Matrix<double, 3, 1>::Unit(2))
        );
        model.TranslateBeam(blade_elem_id, std::array{hub_radius, 0., 0.});
        model.RotateBeamAboutPoint(
            blade_elem_id,
            std::array{
                rotation_quaternion.w(), rotation_quaternion.x(), rotation_quaternion.y(),
                rotation_quaternion.z()
            },
            origin
        );
        model.SetBeamVelocityAboutPoint(blade_elem_id, velocity, origin);
    }

    // Add a prescribed BC for each root node
    for (const auto& beam_element : model.GetBeamElements()) {
        model.AddPrescribedBC(beam_element.node_ids.front());
    }

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(4);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    const double t_end(0.1);
    const auto num_steps = static_cast<size_t>(std::floor(t_end / step_size + 1.0));
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Perform time steps and check for convergence within max_iter iterations
    for (auto i : std::views::iota(0U, num_steps)) {
        // Calculate hub rotation for this time step
        const auto hub_angle = step_size * (i + 1);
        const auto q_hub = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(hub_angle, Eigen::Matrix<double, 3, 1>(&velocity[3]))
        );

        // Define hub translation/rotation displacement
        const auto u_hub = std::array{0., 0., 0., q_hub.w(), q_hub.x(), q_hub.y(), q_hub.z()};

        // Update constraint displacements
        for (auto j : std::views::iota(0U, constraints.num_constraints)) {
            constraints.UpdateDisplacement(j, u_hub);
        }

        // Take step
        auto converged = Step(parameters, solver, elements, state, constraints);

        // Verify that step converged
        EXPECT_EQ(converged, true);
    }
}

TEST(RotatingBeamTest, MasslessConstraints) {
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    constexpr double omega = 0.1;
    constexpr double hub_size = 2.;
    std::vector<size_t> node_ids;
    std::ranges::transform(node_s, std::back_inserter(node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode().SetElemLocation(s).SetPosition(x, 0., 0., 1., 0., 0., 0.).Build();
    });

    // Add beam element and set its position and velocity
    model.AddBeamElement(node_ids, sections, quadrature);
    model.TranslateBeam(0, std::array{hub_size, 0., 0.});
    model.SetBeamVelocityAboutPoint(
        0, std::array{0., 0., 0., 0., 0., omega}, std::array{0., 0., 0.}
    );

    // Add hub node and associated constraints
    auto hub_node_id_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
    model.AddRigidJointConstraint(std::array{hub_node_id_id, node_ids.front()});
    auto hub_bc_id = model.AddPrescribedBC(hub_node_id_id);

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(5);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Perform 10 time steps and check for convergence within max_iter iterations
    for (auto i : std::views::iota(0, 10)) {
        // Set constraint displacement
        const auto hub_angle = omega * step_size * (i + 1);
        const auto q = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(hub_angle, Eigen::Matrix<double, 3, 1>::Unit(2))
        );
        const auto displacement = std::array{0., 0., 0., q.w(), q.x(), q.y(), q.z()};
        constraints.UpdateDisplacement(hub_bc_id, displacement);
        const auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_EQ(converged, true);
    }

    expect_kokkos_view_2D_equal(
        state.q,
        {{-0.000099661884299369481, 0.019999672917628962, -3.6608854058480302E-25,
          0.99998750002604175, -1.5971376141505654E-26, 3.1592454262792375E-25,
          0.004999979166692714},
         {-0.00015838391157346692, 0.031746709275713193, -2.8155520815870626E-13,
          0.99998750002143066, 2.7244338869052949E-12, 1.989181042516661E-12, 0.0049999800888738608},
         {-0.00027859681974392133, 0.055737500699772298, 2.815269319303426E-12, 0.9999875000205457,
          7.3510877107173739E-12, 1.0550370096863904E-12, 0.0049999802658924715},
         {-0.00042131446700509681, 0.08426017738413949, 8.2854411551089936E-12, 0.99998750002161218,
          3.7252296525466957E-11, -5.26890056047209E-14, 0.0049999800525935617},
         {-0.00054093210652801399, 0.10825097509997549, -9.3934322245617647E-12, 0.99998750002142056,
          4.0321076018153484E-11, 5.2579938812420674E-12, 0.0049999800909203019},
         {-0.00059944528351138049, 0.11999801747595988, -2.6207280972097857E-11, 0.99998750002237801,
          3.4435006114567926E-11, 6.4250095159262128E-12, 0.0049999798994432168},
         {0., 0., 0., 0.99998750002604219, 2.2269013449027429E-29, 1.884955233551297E-29,
          0.0049999791666927107}}
    );
}

TEST(RotatingBeamTest, RotationControlConstraint) {
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    constexpr double hub_size = 2.;
    std::vector<size_t> node_ids;
    std::ranges::transform(node_s, std::back_inserter(node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode().SetElemLocation(s).SetPosition(x, 0., 0., 1., 0., 0., 0.).Build();
    });

    // Add beam element and set its position and velocity
    model.AddBeamElement(node_ids, sections, quadrature);
    model.TranslateBeam(0, std::array{hub_size, 0., 0.});

    // Add hub node and associated constraints
    auto pitch = 0.;
    auto hub_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
    model.AddRotationControl(
        std::array{hub_node_id, node_ids.front()}, std::array{1., 0., 0.}, &pitch
    );
    model.AddFixedBC(hub_node_id);

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(5);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Perform 10 time steps and check for convergence within max_iter iterations
    for (auto i : std::views::iota(0, 10)) {
        const auto t = step_size * static_cast<double>(i + 1);
        // Set pitch
        pitch = t * std::numbers::pi / 2.;
        const auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_EQ(converged, true);
    }

    // Check that remaining displacements are as expected
    expect_kokkos_view_2D_equal(
        state.q,
        {{-5.7690945215728628E-18, 2.0652319043893875E-18, -2.4953577261422928E-20,
          0.99691733356165013, 0.078459097906677058, 5.0946025505270351E-19, 3.3980589449407732E-19},
         {-2.9904494403209058E-7, 0.00014453413260242541, -0.00075720167307353353,
          0.99700858667767766, 0.077290133367305239, 0.00033527696532594382,
          0.000031054788779458509},
         {-0.000001465015511032517, 0.00057112956926573543, -0.0031674836974124104,
          0.99711128371780355, 0.075953862519306012, 0.00031282312872539991,
          0.000028155321221201416},
         {-0.0000020860700662326159, 0.00062108133928823466, -0.0035095491673744371,
          0.99720257290191016, 0.074745725742637242, -0.00032312502616037265,
          -0.000025908193827167231},
         {-0.0000043494810453724702, 0.00012641646939721653, -0.00026687208011275229,
          0.99728073405283557, 0.073693469884480778, -0.00063846060120356163,
          -0.000048492096870836455},
         {-0.0000058187577280119813, -0.00014174393267219577, 0.0015608892113182621,
          0.9972883429355086, 0.073590276246437658, -0.00065565396382017294,
          -0.000049038763777556399},
         {0, 0, 0, 1, 0, 0, 0}}
    );
}

TEST(RotatingBeamTest, CompoundRotationControlConstraint) {
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    constexpr double hub_size = 2.;
    std::vector<size_t> node_ids;
    std::ranges::transform(node_s, std::back_inserter(node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode().SetElemLocation(s).SetPosition(x, 0., 0., 1., 0., 0., 0.).Build();
    });

    // Add beam element and set its position and velocity
    model.AddBeamElement(node_ids, sections, quadrature);
    model.TranslateBeam(0, std::array{hub_size, 0., 0.});

    // Add hub node and associated constraints
    auto pitch = 0.;
    auto hub_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
    model.AddRotationControl(std::array{hub_node_id, node_ids[0]}, std::array{1., 0., 0.}, &pitch);
    auto hub_bc_id = model.AddPrescribedBC(hub_node_id);

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(5);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    double azimuth = 0.;

    // Perform 10 time steps and check for convergence within max_iter iterations
    for (auto i : std::views::iota(0, 100)) {
        const auto t = step_size * static_cast<double>(i + 1);
        pitch = t * std::numbers::pi / 2.;
        azimuth = 0.5 * t * std::numbers::pi / 2.;
        const auto q =
            Eigen::Quaternion<double>(Eigen::AngleAxis(azimuth, Eigen::Matrix<double, 3, 1>::Unit(2))
            );
        const auto displacement = std::array{0., 0., 0., q.w(), q.x(), q.y(), q.z()};
        constraints.UpdateDisplacement(hub_bc_id, displacement);
        const auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_EQ(converged, true);
    }

    const auto q = kokkos_view_2D_to_vector(state.q);
    const auto aa =
        Eigen::AngleAxis<double>(Eigen::Quaternion<double>(q[0][3], q[0][4], q[0][5], q[0][6]));
    const auto rv = aa.axis() * aa.angle();

    // Same as euler rotation xz [azimuth, pitch]
    EXPECT_NEAR(rv[0], 1.482189821649821, 1e-8);
    EXPECT_NEAR(rv[1], 0.61394312430788889, 1e-8);
    EXPECT_NEAR(rv[2], 0.61394312416734476, 1e-8);
}

TEST(RotatingBeamTest, RevoluteJointConstraint) {
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    const double omega = 0.1;
    constexpr double hub_size = 2.;
    std::vector<size_t> node_ids;
    std::ranges::transform(node_s, std::back_inserter(node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode().SetElemLocation(s).SetPosition(x, 0., 0., 1., 0., 0., 0.).Build();
    });

    // Add beam element and set its position and velocity
    model.AddBeamElement(node_ids, sections, quadrature);
    model.TranslateBeam(0, std::array{hub_size, 0., 0.});
    model.SetBeamVelocityAboutPoint(
        0, std::array{0., 0., 0., 0., 0., omega}, std::array{0., 0., 0.}
    );

    // Add hub node and ground node
    auto hub_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
    auto ground_node_id = model.AddNode().SetPosition(0, 0., -1., 1., 0., 0., 0.).Build();

    // Add constraints
    model.AddFixedBC(ground_node_id);  // Ground node is fixed

    // Revolute joint constraint
    auto torque = 0.;
    model.AddRevoluteJointConstraint(
        std::array{ground_node_id, hub_node_id}, std::array{0., 0., 0.}, &torque
    );

    // Hub node is rigidly connected
    model.AddRigidJointConstraint(std::array{hub_node_id, node_ids.front()});

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(5);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.9);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Run 10 steps
    for ([[maybe_unused]] auto i : std::views::iota(0, 5)) {
        const auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_EQ(converged, true);
    }

    expect_kokkos_view_2D_equal(
        state.q,
        {{-0.00002499992806604526, 0.009999954363364278, 7.0592758596667977E-26, 0.99999687500410894,
          8.0588517542554213E-29, 5.3983736084468692E-26, 0.0024999964033195574},
         {-0.000039681188471502004, 0.015873544516367064, 4.8213624609892645E-14,
          0.99999687500405054, -1.4117736696803771E-13, 4.0393896390903965E-14,
          0.0024999964266638514},
         {-0.000069668915568583335, 0.027869085208344892, 4.0193427070238621E-14,
          0.99999687500401824, -4.6940123973356305E-13, 7.4011787713126742E-14,
          0.0024999964396126765},
         {-0.00010532153768001662, 0.0421305961297474, 9.1114041847296216E-15, 0.99999687500407363,
          -2.4869679200969056E-13, 6.7014722734775103E-15, 0.0024999964174613749},
         {-0.00013530744593540648, 0.054126136607669233, -5.3259699489140029E-15,
          0.99999687500410672, 1.6008699820757778E-13, -6.0232179820024248E-14,
          0.0024999964042487837},
         {-0.00014999065003689889, 0.059999726696767466, 9.5097860930873942E-14, 0.9999968750041126,
          2.113602318183918E-13, -7.958169516703918E-14, 0.0024999964018790257},
         {0, 0, 0, 0.99999687500410894, 7.2192479817696718E-29, -2.7782381382812961E-27,
          0.0024999964033195574},
         {0, 0, 0, 1, 0, 0, 0}}
    );
}

void GeneratorTorqueWithAxisTilt(
    double tilt, const std::vector<double>& expected_azimuth_q,
    const std::vector<double>& expected_azimuth_vel,
    const std::vector<double>& expected_revolute_joint_output
) {
    auto model = Model();

    // Gravity vector - assume no gravity
    model.SetGravity(0., 0., 0.);

    // Calculate tilt about x axis as a quaternion
    const auto node_tilt = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(tilt, Eigen::Matrix<double, 3, 1>::Unit(0))
    );
    const auto node_tilt_array =
        std::array{node_tilt.w(), node_tilt.x(), node_tilt.y(), node_tilt.z()};

    // Build vector of nodes (straight along x axis, no rotation)
    constexpr double hub_size = 2.;
    std::vector<size_t> node_ids;
    std::ranges::transform(node_s, std::back_inserter(node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode().SetElemLocation(s).SetPosition(x, 0., 0., 1., 0., 0., 0.).Build();
    });

    // Add beam element and set its position and velocity
    model.AddBeamElement(node_ids, sections, quadrature);
    model.TranslateBeam(0, std::array{hub_size, 0., 0.});
    model.RotateBeamAboutPoint(0, node_tilt_array, std::array{0., 0., 0.});

    // Add shaft base, azimuth, and hub nodes as massless points
    auto shaft_base_node_id = model.AddNode().SetPosition(0, 0., 0., 1., 0., 0., 0.).Build();
    auto azimuth_node_id = model.AddNode().SetPosition(0, 0, 0, 1., 0., 0., 0.).Build();
    auto hub_node_id =
        model.AddNode().SetPosition(0, std::sin(tilt), std::cos(tilt), 1., 0., 0., 0.).Build();

    // Add constraints between the nodes to simulate a rotor with a generator
    model.AddFixedBC(shaft_base_node_id);  // Fixed shaft base

    // Add torque to the azimuth node to simulate generator torque
    auto torque = 100.;
    auto shaft_rj_id = model.AddRevoluteJointConstraint(  // Azimuth can rotate around shaft base
        std::array{shaft_base_node_id, azimuth_node_id},
        std::array{0., std::sin(tilt), std::cos(tilt)}, &torque
    );

    // Hub is rigidly attached to azimuth
    model.AddRigidJointConstraint(std::array{azimuth_node_id, hub_node_id});

    // Beam is rigidly attached to hub
    model.AddRigidJointConstraint(std::array{hub_node_id, node_ids.front()});

    // Solution parameters
    const bool is_dynamic_solve(true);
    const int max_iter(5);
    const double step_size(0.01);  // seconds
    const double rho_inf(0.);
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Run 10 steps
    for ([[maybe_unused]] auto i : std::views::iota(0, 10)) {
        const auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_EQ(converged, true);
    }

    // Check that the azimuth node has rotated by the expected amount
    auto azimuth_q = Kokkos::View<double[7]>("azimuth_q");
    Kokkos::deep_copy(azimuth_q, Kokkos::subview(state.q, azimuth_node_id, Kokkos::ALL));
    expect_kokkos_view_1D_equal(azimuth_q, expected_azimuth_q);

    // Check the azimuth node angular velocity is as expected
    auto azimuth_vel = Kokkos::View<double[6]>("azimuth_vel");
    Kokkos::deep_copy(azimuth_vel, Kokkos::subview(state.v, azimuth_node_id, Kokkos::ALL));
    expect_kokkos_view_1D_equal(azimuth_vel, expected_azimuth_vel);

    // Get revolute joint output
    auto revolute_joint_out = Kokkos::View<double[3]>("revolute_joint_out");
    Kokkos::deep_copy(
        revolute_joint_out, Kokkos::subview(constraints.output, shaft_rj_id, Kokkos::ALL)
    );
    // Check output (azimuth, angular velocity, angular acceleration)
    expect_kokkos_view_1D_equal(revolute_joint_out, expected_revolute_joint_output);
}

TEST(RotatingBeamTest, GeneratorTorque_Tilt0) {
    GeneratorTorqueWithAxisTilt(
        0.,                                            // Shaft tilt
        {0., 0., 0., 0.99998634, 0., 0., -0.0052267},  // Azimuth node rotational displacement
        {0., 0., 0., 0., 0., -0.18978539},             // Azimuth node rotational velocity
        {-0.01045353, -0.18978539, -0.566558}          // Shaft angular rotation, velocity, accel
    );
}

TEST(RotatingBeamTest, GeneratorTorque_Tilt90) {
    GeneratorTorqueWithAxisTilt(
        std::numbers::pi / 2.,                         // Shaft tilt
        {0., 0., 0., 0.99998634, 0., -0.0052267, 0.},  // Azimuth node rotational displacement
        {0., 0., 0., 0., -0.18978539, 0.},             // Azimuth node rotational velocity
        {-0.01045353, -0.18978539, -0.566560}          // Shaft angular rotation, velocity, accel
    );
}

}  // namespace kynema::tests
