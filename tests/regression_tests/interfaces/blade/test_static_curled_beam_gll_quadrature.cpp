#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

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

TEST(BladeInterfaceTest, StaticCurledBeam_GllQuadrature_WholeBeam) {
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
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::WholeBeam)
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

TEST(BladeInterfaceTest, StaticCurledBeam_GllQuadrature_ThreeSections_WholeBeam) {
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
        .SetQuadratureStyle(interfaces::components::BeamInput::QuadratureStyle::WholeBeam)
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

    EXPECT_NEAR(tip_positions[4][0], -1.8920340550454853, 1e-8);
    EXPECT_NEAR(tip_positions[4][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[4][2], 1.3745940385447364, 1e-8);

    EXPECT_NEAR(tip_positions[5][0], -2.2205412397724444e-05, 1e-8);
    EXPECT_NEAR(tip_positions[5][1], 0., 1e-8);
    EXPECT_NEAR(tip_positions[5][2], -2.68489535049099e-08, 1e-8);
}
}  // namespace kynema::tests
