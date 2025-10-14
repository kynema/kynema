#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "interfaces/blade/blade_interface.hpp"
#include "interfaces/blade/blade_interface_builder.hpp"
#include "interfaces/components/beam_builder.hpp"

namespace kynema::tests {

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
}  // namespace kynema::tests
