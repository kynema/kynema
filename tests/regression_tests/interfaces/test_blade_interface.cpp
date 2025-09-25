#include <numbers>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"
#include "step/step.hpp"

namespace kynema::tests {

TEST(BladeInterfaceTest, BladeWindIO) {
    // Read WindIO yaml file
    const YAML::Node wio = YAML::LoadFile("interfaces_test_files/IEA-15-240-RWT.yaml");
    const auto& wio_blade = wio["components"]["blade"];
    const auto write_output{true};

    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};

    // Set solution parameters
    const double time_step{0.01};
    builder.Solution()
        .EnableDynamicSolve()
        .SetTimeStep(time_step)
        .SetDampingFactor(0.0)
        .SetMaximumNonlinearIterations(6)
        .SetAbsoluteErrorTolerance(1e-6)
        .SetRelativeErrorTolerance(1e-4);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.BladeWindIO");
    }

    // Set blade parameters
    // Use prescribed root motion to fix root rotation
    const auto n_nodes{11};
    builder.Blade().SetElementOrder(n_nodes - 1).PrescribedRootMotion(true);

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
        builder.Blade().AddRefAxisTwist(twist_grid[i], twist_values[i] * std::numbers::pi / 180.);
    }

    // Add blade section properties
    const auto stiff_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["stiff_matrix"];
    const auto inertia_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["inertia_matrix"];
    const auto k_grid = stiff_matrix["grid"].as<std::vector<double>>();
    const auto m_grid = inertia_matrix["grid"].as<std::vector<double>>();
    const auto n_sections = k_grid.size();
    if (m_grid.size() != k_grid.size()) {
        throw std::runtime_error("stiffness and mass matrices not on same grid");
    }
    for (auto i : std::views::iota(0U, n_sections)) {
        if (abs(m_grid[i] - k_grid[i]) > 1e-8) {
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

    // Build blade interface
    auto interface = builder.Build();

    // Get reference to tip node
    auto& tip_node = interface.Blade().nodes.back();

    // Loop through solution iterations
    for (auto i : std::views::iota(1U, 20U)) {
        // Calculate time
        const auto t{static_cast<double>(i) * time_step};

        // Apply oscillating moment about z axis to tip
        tip_node.loads[1] = 2.0e5 * sin(2. * std::numbers::pi * 0.05 * t);

        // Take step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);
    }

    EXPECT_NEAR(tip_node.position[0], 117.28591620612008, 1e-10);
    EXPECT_NEAR(tip_node.position[1], 0.1714518799428682, 1e-10);
    EXPECT_NEAR(tip_node.position[2], 4.0011349824240705, 1e-10);
    EXPECT_NEAR(tip_node.position[3], 0.99876122551364266, 1e-10);
    EXPECT_NEAR(tip_node.position[4], -0.003559701950531693, 1e-10);
    EXPECT_NEAR(tip_node.position[5], -0.049343778710499483, 1e-10);
    EXPECT_NEAR(tip_node.position[6], 0.0053417632931122526, 1e-10);
}

TEST(BladeInterfaceTest, RotatingBeam) {
    const auto time_step{0.01};
    const auto omega = std::array{0., 0., 1.};
    const auto x0_root = std::array{2., 0., 0.};
    const auto root_vel = math::CrossProduct(omega, x0_root);
    const auto write_output{false};

    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};

    // Set solution options
    builder.Solution()
        .EnableDynamicSolve()
        .SetTimeStep(time_step)
        .SetDampingFactor(0.0)
        .SetMaximumNonlinearIterations(6)
        .SetAbsoluteErrorTolerance(1e-6)
        .SetRelativeErrorTolerance(1e-4);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.RotatingBeam");
    }

    // Node locations (GLL quadrature)
    const auto node_s = std::vector{
        0., 0.11747233803526763, 0.35738424175967748, 0.64261575824032247, 0.88252766196473242, 1.
    };

    builder.Blade()
        .SetElementOrder(5)
        .PrescribedRootMotion(true)
        .SetRootPosition({x0_root[0], x0_root[1], x0_root[2], 1., 0., 0., 0.})
        .SetRootVelocity({root_vel[0], root_vel[1], root_vel[2], omega[0], omega[1], omega[2]})
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    // Add reference axis coordinates and twist
    for (const auto s : node_s) {
        builder.Blade().AddRefAxisPoint(
            s, {10. * s, 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{
        0.0254460438286208, 0.1292344072003028, 0.2970774243113015, 0.5,
        0.7029225756886985, 0.8707655927996972, 0.9745539561713792,
    };

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            std::array{
                std::array{8.538e-2, 0., 0., 0., 0., 0.},
                std::array{0., 8.538e-2, 0., 0., 0., 0.},
                std::array{0., 0., 8.538e-2, 0., 0., 0.},
                std::array{0., 0., 0., 1.4433e-2, 0., 0.},
                std::array{0., 0., 0., 0., 0.40972e-2, 0.},
                std::array{0., 0., 0., 0., 0., 1.0336e-2},
            },
            std::array{
                std::array{1368.17e3, 0., 0., 0., 0., 0.},
                std::array{0., 88.56e3, 0., 0., 0., 0.},
                std::array{0., 0., 38.78e3, 0., 0., 0.},
                std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
                std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
                std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Create the interface
    auto interface = builder.Build();

    // Loop through time steps
    for (auto i : std::views::iota(1U, 10U)) {
        // Calculate time at end of step
        const auto t{static_cast<double>(i) * time_step};

        // Calculate root displacement from initial position and apply
        const auto u_rot =
            math::RotationVectorToQuaternion({omega[0] * t, omega[1] * t, omega[2] * t});
        const auto x_root = math::RotateVectorByQuaternion(u_rot, x0_root);
        const auto u_trans =
            std::array{x_root[0] - x0_root[0], x_root[1] - x0_root[1], x_root[2] - x0_root[2]};
        interface.SetRootDisplacement(
            {u_trans[0], u_trans[1], u_trans[2], u_rot[0], u_rot[1], u_rot[2], u_rot[3]}
        );

        // Calculate state at end of step
        const auto converged = interface.Step();

        // Check for convergence
        ASSERT_EQ(converged, true);
    }

    //--------------------------------------------------------------------------
    // Root Node
    //--------------------------------------------------------------------------

    EXPECT_NEAR(interface.Blade().nodes[0].position[0], 1.9919054660239885, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].position[1], 0.17975709839602211, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].position[2], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].position[3], 0.99898767084784234, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].position[4], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].position[5], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].position[6], 0.044984814037660227, 1e-10);

    EXPECT_NEAR(interface.Blade().nodes[0].displacement[0], -0.0080945339760114532, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].displacement[1], 0.17975709839602211, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].displacement[2], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].displacement[3], 0.99898767084784234, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].displacement[4], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].displacement[5], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].displacement[6], 0.044984814037660227, 1e-10);

    EXPECT_NEAR(interface.Blade().nodes[0].velocity[0], -0.17976259212149071, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].velocity[1], 1.9919719054875213, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].velocity[2], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].velocity[3], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].velocity[4], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].velocity[5], 1., 1e-10);

    EXPECT_NEAR(interface.Blade().nodes[0].acceleration[0], -1.9920882238788629, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].acceleration[1], -0.17977158309645525, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].acceleration[2], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].acceleration[3], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].acceleration[4], 0., 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[0].acceleration[5], 0., 1e-10);

    //--------------------------------------------------------------------------
    // Tip Node
    //--------------------------------------------------------------------------

    EXPECT_NEAR(interface.Blade().nodes[5].position[0], 11.951460424343411, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].position[1], 1.0785463019560009, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].position[2], 2.2845605555050727E-10, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].position[3], 0.99898766738195299, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].position[4], 2.6570959740959799E-9, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].position[5], -1.0123941616837804E-10, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].position[6], 0.044984891005365443, 1e-10);

    EXPECT_NEAR(interface.Blade().nodes[5].displacement[0], -0.048539575656588806, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].displacement[1], 1.0785463019560009, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].displacement[2], 2.2845605555050727E-10, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].displacement[3], 0.99898766738195299, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].displacement[4], 2.6570959740959799E-9, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].displacement[5], -1.0123941616837804E-10, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].displacement[6], 0.044984891005365443, 1e-10);

    EXPECT_NEAR(interface.Blade().nodes[5].velocity[0], -1.0786835619583102, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].velocity[1], 11.95303967174639, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].velocity[2], -6.3862048021175104E-8, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].velocity[3], 0.000003159329892760129, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].velocity[4], 2.5866107917084571E-7, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].velocity[5], 1.0001562024381443, 1e-10);

    EXPECT_NEAR(interface.Blade().nodes[5].acceleration[0], -11.955952417551353, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].acceleration[1], -1.0800839943831992, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].acceleration[2], -0.000013544534439173783, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].acceleration[3], -0.000046709402446003812, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].acceleration[4], 0.000008256203749646939, 1e-10);
    EXPECT_NEAR(interface.Blade().nodes[5].acceleration[5], 0.00020635762959181575, 1e-10);
}

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

TEST(BladeInterfaceTest, TwoBeams) {
    // Create interface builder
    auto builder = interfaces::components::BeamBuilder{};

    const auto n_nodes = 4;

    builder.SetElementOrder(n_nodes - 1)
        .PrescribedRootMotion(false)
        .SetRootPosition({0., 0., 0., 1., 0., 0., 0.})
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    // Add reference axis coordinates and twist
    const auto n_kps = 12;
    for (auto i : std::views::iota(0, n_kps)) {
        const auto s = static_cast<double>(i) / (n_kps - 1);
        builder.AddRefAxisPoint(
            s, {10. * s, 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{0.,   0.05, 0.1,  0.15, 0.2,  0.25, 0.3,
                                        0.35, 0.4,  0.45, 0.5,  0.55, 0.6,  0.65,
                                        0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.};

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.AddSection(
            s,
            std::array{
                std::array{8.538e-2, 0., 0., 0., 0., 0.},
                std::array{0., 8.538e-2, 0., 0., 0., 0.},
                std::array{0., 0., 8.538e-2, 0., 0., 0.},
                std::array{0., 0., 0., 1.4433e-2, 0., 0.},
                std::array{0., 0., 0., 0., 0.40972e-2, 0.},
                std::array{0., 0., 0., 0., 0., 1.0336e-2},
            },
            std::array{
                std::array{1368.17e3, 0., 0., 0., 0., 0.},
                std::array{0., 88.56e3, 0., 0., 0., 0.},
                std::array{0., 0., 38.78e3, 0., 0., 0.},
                std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
                std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
                std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    Model model({0., 0., -9.81});

    auto beam_1 = builder.Build(model);
    auto beam_2 = builder.Build(model);

    for (auto& node : beam_2.nodes) {
        model.GetNode(node.id).x0[0] += 10.;
    }

    model.AddPrescribedBC(beam_1.nodes[0].id);
    model.AddRigidJointConstraint({beam_1.nodes[n_nodes - 1].id, beam_2.nodes[0].id});

    // Create solver parameters
    auto parameters = StepParameters(true, 6, 0.01, 0.);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    for ([[maybe_unused]] auto iteration : std::views::iota(0, 1000)) {
        const auto converged = Step(parameters, solver, elements, state, constraints);
        ASSERT_TRUE(converged);
    }
}

TEST(BladeInterfaceTest, StaticCurledBeam_TrapezoidalQuadrature) {
    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    builder.Solution()
        .EnableStaticSolve()
        .SetTimeStep(1.)
        .SetDampingFactor(1.)
        .SetMaximumNonlinearIterations(10)
        .SetAbsoluteErrorTolerance(1e-5)
        .SetRelativeErrorTolerance(1e-3);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.StaticCurledBeam_Trap");
    }

    // Node locations
    const std::vector<double> kp_s{0., 1.};

    builder.Blade()
        .SetElementOrder(10)
        .SetSectionRefinement(0)
        .PrescribedRootMotion(true)
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{0.,  .05, 0.1, .15, 0.2, .25, 0.3, .35, 0.4, .45, 0.5,
                                        .55, 0.6, .65, 0.7, .75, 0.8, .85, 0.9, .95, 1.};

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            std::array{
                std::array{1., 0., 0., 0., 0., 0.},
                std::array{0., 1., 0., 0., 0., 0.},
                std::array{0., 0., 1., 0., 0., 0.},
                std::array{0., 0., 0., 1., 0., 0.},
                std::array{0., 0., 0., 0., 1., 0.},
                std::array{0., 0., 0., 0., 0., 1.},
            },
            std::array{
                std::array{1770.e3, 0., 0., 0., 0., 0.},
                std::array{0., 1770.e3, 0., 0., 0., 0.},
                std::array{0., 0., 1770.e3, 0., 0., 0.},
                std::array{0., 0., 0., 8.16e3, 0., 0.},
                std::array{0., 0., 0., 0., 86.9e3, 0.},
                std::array{0., 0., 0., 0., 0., 215.e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    // Create vector to store deformed tip positions
    std::vector<std::array<double, 3>> tip_positions;

    // Get reference to tip node
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];

    // Loop through moments to apply to tip
    const std::vector<double> moments{0., 10920.0, 21840.0, 32761.0, 43681.0, 54601.0};
    for (const auto m : moments) {
        // Apply moment to tip about y axis
        tip_node.loads[4] = -m;

        // Static step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);

        // Add tip position
        tip_positions.emplace_back(
            std::array{tip_node.position[0], tip_node.position[1], tip_node.position[2]}
        );
    }

    EXPECT_NEAR(tip_positions[0][0], 10., 1e-8);
    EXPECT_NEAR(tip_positions[0][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[0][2], 0., 1e-8);

    EXPECT_NEAR(tip_positions[1][0], 7.5396813678794645, 1e-8);
    EXPECT_NEAR(tip_positions[1][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[1][2], 5.5363879677563252, 1e-8);

    EXPECT_NEAR(tip_positions[2][0], 2.275087482106132, 1e-8);
    EXPECT_NEAR(tip_positions[2][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[2][2], 7.2166560706481686, 1e-8);

    EXPECT_NEAR(tip_positions[3][0], -1.6222827675944771, 1e-8);
    EXPECT_NEAR(tip_positions[3][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[3][2], 4.7694966546892239, 1e-8);

    EXPECT_NEAR(tip_positions[4][0], -1.9054629771623546, 1e-8);
    EXPECT_NEAR(tip_positions[4][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[4][2], 1.3243726828814482, 1e-8);

    EXPECT_NEAR(tip_positions[5][0], 0.021386893541979646, 1e-8);
    EXPECT_NEAR(tip_positions[5][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[5][2], 0.0006097054603659835, 1e-8);
}

TEST(BladeInterfaceTest, StaticCurledBeam_GllQuadrature) {
    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    builder.Solution()
        .EnableStaticSolve()
        .SetTimeStep(1.)
        .SetDampingFactor(1.)
        .SetMaximumNonlinearIterations(10)
        .SetAbsoluteErrorTolerance(1e-5)
        .SetRelativeErrorTolerance(1e-3);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.StaticCurledBeam_GLL");
    }

    // Node locations
    const std::vector<double> kp_s{0., 1.};

    builder.Blade()
        .SetElementOrder(10)
        .SetSectionRefinement(14)
        .PrescribedRootMotion(true)
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{0., 1.};

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            std::array{
                std::array{1., 0., 0., 0., 0., 0.},
                std::array{0., 1., 0., 0., 0., 0.},
                std::array{0., 0., 1., 0., 0., 0.},
                std::array{0., 0., 0., 1., 0., 0.},
                std::array{0., 0., 0., 0., 1., 0.},
                std::array{0., 0., 0., 0., 0., 1.},
            },
            std::array{
                std::array{1770.e3, 0., 0., 0., 0., 0.},
                std::array{0., 1770.e3, 0., 0., 0., 0.},
                std::array{0., 0., 1770.e3, 0., 0., 0.},
                std::array{0., 0., 0., 8.16e3, 0., 0.},
                std::array{0., 0., 0., 0., 86.9e3, 0.},
                std::array{0., 0., 0., 0., 0., 215.e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    // Create vector to store deformed tip positions
    std::vector<std::array<double, 3>> tip_positions;

    // Get reference to tip node
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];

    // Loop through moments to apply to tip
    const std::vector<double> moments{0., 10920.0, 21840.0, 32761.0, 43681.0, 54601.0};
    for (const auto m : moments) {
        // Apply moment to tip about y axis
        tip_node.loads[4] = -m;

        // Static step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);

        // Add tip position
        tip_positions.emplace_back(
            std::array{tip_node.position[0], tip_node.position[1], tip_node.position[2]}
        );
    }

    EXPECT_NEAR(tip_positions[0][0], 10., 1e-8);
    EXPECT_NEAR(tip_positions[0][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[0][2], 0., 1e-8);

    EXPECT_NEAR(tip_positions[1][0], 7.568339485511153, 1e-8);
    EXPECT_NEAR(tip_positions[1][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[1][2], 5.4986033633097984, 1e-8);

    EXPECT_NEAR(tip_positions[2][0], 2.3388913532848772, 1e-8);
    EXPECT_NEAR(tip_positions[2][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[2][2], 7.1978711995719049, 1e-8);

    EXPECT_NEAR(tip_positions[3][0], -1.5592428635408258, 1e-8);
    EXPECT_NEAR(tip_positions[3][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[3][2], 4.7984128627696361, 1e-8);

    EXPECT_NEAR(tip_positions[4][0], -1.8920340550454853, 1e-8);
    EXPECT_NEAR(tip_positions[4][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[4][2], 1.3745940385447364, 1e-8);

    EXPECT_NEAR(tip_positions[5][0], -2.2205412397724444e-05, 1e-8);
    EXPECT_NEAR(tip_positions[5][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[5][2], -2.68489535049099e-08, 1e-8);
}

TEST(BladeInterfaceTest, StaticCurledBeam_GllQuadrature_ThreeSections) {
    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    builder.Solution()
        .EnableStaticSolve()
        .SetTimeStep(1.)
        .SetDampingFactor(1.)
        .SetMaximumNonlinearIterations(10)
        .SetAbsoluteErrorTolerance(1e-5)
        .SetRelativeErrorTolerance(1e-3);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.StaticCurledBeam_GLL");
    }

    // Node locations
    const std::vector<double> kp_s{0., 1.};

    builder.Blade()
        .SetElementOrder(10)
        .SetSectionRefinement(7)
        .PrescribedRootMotion(true)
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{0., 0.5, 1.};

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            std::array{
                std::array{1., 0., 0., 0., 0., 0.},
                std::array{0., 1., 0., 0., 0., 0.},
                std::array{0., 0., 1., 0., 0., 0.},
                std::array{0., 0., 0., 1., 0., 0.},
                std::array{0., 0., 0., 0., 1., 0.},
                std::array{0., 0., 0., 0., 0., 1.},
            },
            std::array{
                std::array{1770.e3, 0., 0., 0., 0., 0.},
                std::array{0., 1770.e3, 0., 0., 0., 0.},
                std::array{0., 0., 1770.e3, 0., 0., 0.},
                std::array{0., 0., 0., 8.16e3, 0., 0.},
                std::array{0., 0., 0., 0., 86.9e3, 0.},
                std::array{0., 0., 0., 0., 0., 215.e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    // Create vector to store deformed tip positions
    std::vector<std::array<double, 3>> tip_positions;

    // Get reference to tip node
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];

    // Loop through moments to apply to tip
    const std::vector<double> moments{0., 10920.0, 21840.0, 32761.0, 43681.0, 54601.0};
    for (const auto m : moments) {
        // Apply moment to tip about y axis
        tip_node.loads[4] = -m;

        // Static step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);

        // Add tip position
        tip_positions.emplace_back(
            std::array{tip_node.position[0], tip_node.position[1], tip_node.position[2]}
        );
    }

    EXPECT_NEAR(tip_positions[0][0], 10., 1e-8);
    EXPECT_NEAR(tip_positions[0][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[0][2], 0., 1e-8);

    EXPECT_NEAR(tip_positions[1][0], 7.568339485511153, 1e-8);
    EXPECT_NEAR(tip_positions[1][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[1][2], 5.4986033633097984, 1e-8);

    EXPECT_NEAR(tip_positions[2][0], 2.3388913532848772, 1e-8);
    EXPECT_NEAR(tip_positions[2][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[2][2], 7.1978711995719049, 1e-8);

    EXPECT_NEAR(tip_positions[3][0], -1.5592428635408258, 1e-8);
    EXPECT_NEAR(tip_positions[3][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[3][2], 4.7984128627696361, 1e-8);

    EXPECT_NEAR(tip_positions[4][0], -1.8920341370666414, 1e-8);
    EXPECT_NEAR(tip_positions[4][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[4][2], 1.3745940986394585, 1e-8);

    EXPECT_NEAR(tip_positions[5][0], -2.6492263431876495e-05, 1e-8);
    EXPECT_NEAR(tip_positions[5][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[5][2], -1.1189941773181999e-07, 1e-8);
}

TEST(BladeInterfaceTest, StaticCurledBeam_GlQuadrature) {
    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    builder.Solution()
        .EnableStaticSolve()
        .SetTimeStep(1.)
        .SetDampingFactor(1.)
        .SetMaximumNonlinearIterations(10)
        .SetAbsoluteErrorTolerance(1e-5)
        .SetRelativeErrorTolerance(1e-3);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.StaticCurledBeam_GL");
    }

    // Node locations
    const std::vector<double> kp_s{0., 1.};

    builder.Blade()
        .SetElementOrder(10)
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetSectionRefinement(14)
        .PrescribedRootMotion(true)
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{0., 1.};

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            std::array{
                std::array{1., 0., 0., 0., 0., 0.},
                std::array{0., 1., 0., 0., 0., 0.},
                std::array{0., 0., 1., 0., 0., 0.},
                std::array{0., 0., 0., 1., 0., 0.},
                std::array{0., 0., 0., 0., 1., 0.},
                std::array{0., 0., 0., 0., 0., 1.},
            },
            std::array{
                std::array{1770.e3, 0., 0., 0., 0., 0.},
                std::array{0., 1770.e3, 0., 0., 0., 0.},
                std::array{0., 0., 1770.e3, 0., 0., 0.},
                std::array{0., 0., 0., 8.16e3, 0., 0.},
                std::array{0., 0., 0., 0., 86.9e3, 0.},
                std::array{0., 0., 0., 0., 0., 215.e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    // Create vector to store deformed tip positions
    std::vector<std::array<double, 3>> tip_positions;

    // Get reference to tip node
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];

    // Loop through moments to apply to tip
    const std::vector<double> moments{0., 10920.0, 21840.0, 32761.0, 43681.0, 54601.0};
    for (const auto m : moments) {
        // Apply moment to tip about y axis
        tip_node.loads[4] = -m;

        // Static step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);

        // Add tip position
        tip_positions.emplace_back(
            std::array{tip_node.position[0], tip_node.position[1], tip_node.position[2]}
        );
    }

    EXPECT_NEAR(tip_positions[0][0], 10., 1e-8);
    EXPECT_NEAR(tip_positions[0][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[0][2], 0., 1e-8);

    EXPECT_NEAR(tip_positions[1][0], 7.568339485511153, 1e-8);
    EXPECT_NEAR(tip_positions[1][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[1][2], 5.4986033633097984, 1e-8);

    EXPECT_NEAR(tip_positions[2][0], 2.3388913532848772, 1e-8);
    EXPECT_NEAR(tip_positions[2][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[2][2], 7.1978711995719049, 1e-8);

    EXPECT_NEAR(tip_positions[3][0], -1.5592428635408258, 1e-8);
    EXPECT_NEAR(tip_positions[3][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[3][2], 4.7984128627696361, 1e-8);

    EXPECT_NEAR(tip_positions[4][0], -1.8920340550454853, 1e-8);
    EXPECT_NEAR(tip_positions[4][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[4][2], 1.3745940385447364, 1e-8);

    EXPECT_NEAR(tip_positions[5][0], -2.2205412397724444e-05, 1e-8);
    EXPECT_NEAR(tip_positions[5][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[5][2], -2.68489535049099e-08, 1e-8);
}

TEST(BladeInterfaceTest, StaticCurledBeam_GlQuadrature_ThreeSections) {
    // Create interface builder
    auto builder = interfaces::BladeInterfaceBuilder{};
    const auto write_output{false};

    builder.Solution()
        .EnableStaticSolve()
        .SetTimeStep(1.)
        .SetDampingFactor(1.)
        .SetMaximumNonlinearIterations(10)
        .SetAbsoluteErrorTolerance(1e-5)
        .SetRelativeErrorTolerance(1e-3);

    if (write_output) {
        builder.Solution().SetOutputFile("BladeInterfaceTest.StaticCurledBeam_GLL");
    }

    // Node locations
    const std::vector<double> kp_s{0., 1.};

    builder.Blade()
        .SetElementOrder(10)
        .SetQuadratureRule(interfaces::components::BeamInput::QuadratureRule::GaussLegendre)
        .SetSectionRefinement(7)
        .PrescribedRootMotion(true)
        .AddRefAxisTwist(0., 0.)
        .AddRefAxisTwist(1., 0.);

    for (const auto s : kp_s) {
        builder.Blade().AddRefAxisPoint(
            s, {s * 10., 0., 0.}, interfaces::components::ReferenceAxisOrientation::X
        );
    }

    // Beam section locations
    const std::vector<double> section_s{0., 0.5, 1.};

    // Add reference axis coordinates and twist
    for (const auto s : section_s) {
        builder.Blade().AddSection(
            s,
            std::array{
                std::array{1., 0., 0., 0., 0., 0.},
                std::array{0., 1., 0., 0., 0., 0.},
                std::array{0., 0., 1., 0., 0., 0.},
                std::array{0., 0., 0., 1., 0., 0.},
                std::array{0., 0., 0., 0., 1., 0.},
                std::array{0., 0., 0., 0., 0., 1.},
            },
            std::array{
                std::array{1770.e3, 0., 0., 0., 0., 0.},
                std::array{0., 1770.e3, 0., 0., 0., 0.},
                std::array{0., 0., 1770.e3, 0., 0., 0.},
                std::array{0., 0., 0., 8.16e3, 0., 0.},
                std::array{0., 0., 0., 0., 86.9e3, 0.},
                std::array{0., 0., 0., 0., 0., 215.e3},
            },
            interfaces::components::ReferenceAxisOrientation::X
        );
    }

    auto interface = builder.Build();

    // Create vector to store deformed tip positions
    std::vector<std::array<double, 3>> tip_positions;

    // Get reference to tip node
    auto& tip_node = interface.Blade().nodes[interface.Blade().nodes.size() - 1];

    // Loop through moments to apply to tip
    const std::vector<double> moments{0., 10920.0, 21840.0, 32761.0, 43681.0, 54601.0};
    for (const auto m : moments) {
        // Apply moment to tip about y axis
        tip_node.loads[4] = -m;

        // Static step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);

        // Add tip position
        tip_positions.emplace_back(
            std::array{tip_node.position[0], tip_node.position[1], tip_node.position[2]}
        );
    }

    EXPECT_NEAR(tip_positions[0][0], 10., 1e-8);
    EXPECT_NEAR(tip_positions[0][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[0][2], 0., 1e-8);

    EXPECT_NEAR(tip_positions[1][0], 7.568339485511153, 1e-8);
    EXPECT_NEAR(tip_positions[1][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[1][2], 5.4986033633097984, 1e-8);

    EXPECT_NEAR(tip_positions[2][0], 2.3388913532848772, 1e-8);
    EXPECT_NEAR(tip_positions[2][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[2][2], 7.1978711995719049, 1e-8);

    EXPECT_NEAR(tip_positions[3][0], -1.5592428635408258, 1e-8);
    EXPECT_NEAR(tip_positions[3][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[3][2], 4.7984128627696361, 1e-8);

    EXPECT_NEAR(tip_positions[4][0], -1.892033981830739, 1e-8);
    EXPECT_NEAR(tip_positions[4][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[4][2], 1.3745939851717117, 1e-8);

    EXPECT_NEAR(tip_positions[5][0], -1.8202489757257467e-05, 1e-8);
    EXPECT_NEAR(tip_positions[5][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[5][2], 6.2366226405430325e-08, 1e-8);
}
}  // namespace kynema::tests
