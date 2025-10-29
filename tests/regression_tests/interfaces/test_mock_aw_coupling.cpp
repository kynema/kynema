#include <array>
#include <numbers>
#include <ranges>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "elements/beams/hollow_circle_properties.hpp"
#include "interfaces/components/aerodynamics.hpp"
#include "interfaces/components/inflow.hpp"
#include "interfaces/turbine/turbine_interface.hpp"
#include "interfaces/turbine/turbine_interface_builder.hpp"

namespace kynema::tests {

TEST(TurbineInterfaceTest, MockAW) {
    const auto time_step{0.005};      // Time step for the simulation
    const auto n_blades{3U};         // Number of blades in turbine
    const auto n_blade_nodes{11};    // Number of nodes per blade
    const auto n_tower_nodes{11};    // Number of nodes in tower
    const auto write_output{true};  // Write output file
    const auto n_steps{100U};

    constexpr auto fluid_density{1.225};
    constexpr auto v_inflow = std::array{6., 0., 0.};

    // Create interface builder
    auto builder = interfaces::TurbineInterfaceBuilder{};

    // Set solution parameters
    builder.Solution()
        .EnableDynamicSolve()
        .SetTimeStep(time_step)
        .SetDampingFactor(0.0)
        .SetGravity({0., 0., -9.81})
        .SetMaximumNonlinearIterations(12)
        .SetAbsoluteErrorTolerance(1e-7)
        .SetRelativeErrorTolerance(1e-6);

    if (write_output) {
        builder.Solution().SetOutputFile("TurbineInterfaceTest.MockAW");
    }

    // Read WindIO yaml file
    const YAML::Node wio = YAML::LoadFile("interfaces_test_files/IEA-15-240-RWT-aero50.yaml");

    // WindIO components
    const auto& wio_tower = wio["components"]["tower"];
    const auto& wio_nacelle = wio["components"]["nacelle"];
    const auto& wio_blade = wio["components"]["blade"];
    const auto& wio_hub = wio["components"]["hub"];

    //--------------------------------------------------------------------------
    // Build Turbine
    //--------------------------------------------------------------------------

    // Get turbine builder
    auto& turbine_builder = builder.Turbine();
    turbine_builder.SetAzimuthAngle(0.)
        .SetRotorApexToHub(0.)
        .SetHubDiameter(wio_hub["diameter"].as<double>())
        .SetConeAngle(wio_hub["cone_angle"].as<double>())
        .SetShaftTiltAngle(wio_nacelle["drivetrain"]["uptilt"].as<double>())
        .SetTowerAxisToRotorApex(wio_nacelle["drivetrain"]["overhang"].as<double>())
        .SetTowerTopToRotorApex(wio_nacelle["drivetrain"]["distance_tt_hub"].as<double>());

    //--------------------------------------------------------------------------
    // Build Blades
    //--------------------------------------------------------------------------

    // Loop through blades and set parameters
    for (auto j : std::views::iota(0U, n_blades)) {
        // Get the blade builder
        auto& blade_builder = turbine_builder.Blade(j);

        // Set blade parameters
        blade_builder.SetElementOrder(n_blade_nodes - 1).PrescribedRootMotion(false);

        // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
        const auto ref_axis = wio_blade["outer_shape_bem"]["reference_axis"];
        const auto axis_grid = ref_axis["x"]["grid"].as<std::vector<double>>();
        const auto x_values = ref_axis["x"]["values"].as<std::vector<double>>();
        const auto y_values = ref_axis["y"]["values"].as<std::vector<double>>();
        const auto z_values = ref_axis["z"]["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, axis_grid.size())) {
            blade_builder.AddRefAxisPoint(
                axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                interfaces::components::ReferenceAxisOrientation::Z
            );
        }

        // Add reference axis twist
        const auto twist = wio_blade["outer_shape_bem"]["twist"];
        const auto twist_grid = twist["grid"].as<std::vector<double>>();
        const auto twist_values = twist["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, twist_grid.size())) {
            blade_builder.AddRefAxisTwist(twist_grid[i], twist_values[i] * std::numbers::pi / 180.);
        }

        // Add blade section properties
        const auto stiff_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["stiff_matrix"];
        const auto inertia_matrix =
            wio_blade["elastic_properties_mb"]["six_x_six"]["inertia_matrix"];
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
            blade_builder.AddSection(
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
    }

    //--------------------------------------------------------------------------
    // Build Tower
    //--------------------------------------------------------------------------

    // Get the tower builder
    auto& tower_builder = turbine_builder.Tower();

    // Set tower parameters
    tower_builder
        .SetElementOrder(n_tower_nodes - 1)  // Set element order to num nodes -1
        .PrescribedRootMotion(false);        // Fix displacement of tower base node

    // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
    const auto t_ref_axis = wio_tower["outer_shape_bem"]["reference_axis"];
    const auto axis_grid = t_ref_axis["x"]["grid"].as<std::vector<double>>();
    const auto x_values = t_ref_axis["x"]["values"].as<std::vector<double>>();
    const auto y_values = t_ref_axis["y"]["values"].as<std::vector<double>>();
    const auto z_values = t_ref_axis["z"]["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, axis_grid.size())) {
        tower_builder.AddRefAxisPoint(
            axis_grid[i], {x_values[i], y_values[i], z_values[i]},
            interfaces::components::ReferenceAxisOrientation::Z
        );
    }

    // Set tower base position from first reference axis point
    const auto tower_base_position =
        std::array<double, 7>{x_values[0], y_values[0], z_values[0], 1., 0., 0., 0.};
    turbine_builder.SetTowerBasePosition(tower_base_position);

    // Add reference axis twist (zero for tower)
    tower_builder.AddRefAxisTwist(0.0, 0.0).AddRefAxisTwist(1.0, 0.0);

    // Find the tower material properties
    const auto t_layer = wio_tower["internal_structure_2d_fem"]["layers"][0];
    const auto t_material_name = t_layer["material"].as<std::string>();
    YAML::Node t_material;
    for (const auto& m : wio["materials"]) {
        if (m["name"] && m["name"].as<std::string>() == t_material_name) {
            t_material = m.as<YAML::Node>();
            break;
        }
    }
    if (!t_material) {
        throw std::runtime_error(
            "Material '" + t_material_name + "' not found in materials section"
        );
    }

    // Add tower section properties
    const auto t_diameter = wio_tower["outer_shape_bem"]["outer_diameter"];
    const auto t_diameter_grid = t_diameter["grid"].as<std::vector<double>>();
    const auto t_diameter_values = t_diameter["values"].as<std::vector<double>>();
    const auto t_wall_thickness = t_layer["thickness"]["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, t_diameter_grid.size())) {
        // Create section mass and stiffness matrices
        const auto section = beams::GenerateHollowCircleSection(
            t_diameter_grid[i], t_material["E"].as<double>(), t_material["G"].as<double>(),
            t_material["rho"].as<double>(), t_diameter_values[i], t_wall_thickness[i],
            t_material["nu"].as<double>()
        );

        // Add section
        tower_builder.AddSection(
            t_diameter_grid[i], section.M_star, section.C_star,
            interfaces::components::ReferenceAxisOrientation::Z
        );
    }

    //--------------------------------------------------------------------------
    // Add mass elements
    //--------------------------------------------------------------------------

    // Get nacelle mass properties from WindIO
    const auto& nacelle_props = wio_nacelle["elastic_properties_mb"];
    const auto system_mass = nacelle_props["system_mass"].as<double>();
    const auto yaw_mass = nacelle_props["yaw_mass"].as<double>();
    const auto system_inertia_tt = nacelle_props["system_inertia_tt"].as<std::vector<double>>();

    // Construct 6x6 inertia matrix for yaw bearing node
    const auto total_mass = system_mass + yaw_mass;
    const auto nacelle_inertia_matrix = std::array<std::array<double, 6>, 6>{
        {{total_mass, 0., 0., 0., 0., 0.},
         {0., total_mass, 0., 0., 0., 0.},
         {0., 0., total_mass, 0., 0., 0.},
         {0., 0., 0., system_inertia_tt[0], system_inertia_tt[3], system_inertia_tt[4]},
         {0., 0., 0., system_inertia_tt[3], system_inertia_tt[1], system_inertia_tt[5]},
         {0., 0., 0., system_inertia_tt[4], system_inertia_tt[5], system_inertia_tt[2]}}
    };

    // Set the nacelle inertia matrix in the turbine builder
    turbine_builder.SetYawBearingInertiaMatrix(nacelle_inertia_matrix);

    // Get hub mass properties from WindIO
    const auto& hub_props = wio_hub["elastic_properties_mb"];
    const auto hub_mass = hub_props["system_mass"].as<double>();
    const auto hub_inertia = hub_props["system_inertia"].as<std::vector<double>>();

    // Construct 6x6 inertia matrix for hub node
    const auto hub_inertia_matrix = std::array<std::array<double, 6>, 6>{
        {{hub_mass, 0., 0., 0., 0., 0.},
         {0., hub_mass, 0., 0., 0., 0.},
         {0., 0., hub_mass, 0., 0., 0.},
         {0., 0., 0., hub_inertia[0], hub_inertia[3], hub_inertia[4]},
         {0., 0., 0., hub_inertia[3], hub_inertia[1], hub_inertia[5]},
         {0., 0., 0., hub_inertia[4], hub_inertia[5], hub_inertia[2]}}
    };

    // Set the hub inertia matrix in the turbine builder
    turbine_builder.SetHubInertiaMatrix(hub_inertia_matrix);

    //--------------------------------------------------------------------------
    // Build Aerodynamics
    //--------------------------------------------------------------------------

    auto& aero_builder =
        builder.Aerodynamics().EnableAero().SetNumberOfAirfoils(1UL).SetAirfoilToBladeMap(
            std::array{0UL, 0UL, 0UL}
        );

    {
        const auto& airfoil_io = wio["airfoils"];
        auto aero_sections = std::vector<interfaces::components::AerodynamicSection>{};
        auto id = 0UL;
        for (const auto& af : airfoil_io) {
            const auto s = af["spanwise_position"].as<double>();
            const auto chord = af["chord"].as<double>();
            const auto twist = af["twist"].as<double>() * std::numbers::pi / 180.;
            const auto section_offset_x = af["section_offset_x"].as<double>();
            const auto section_offset_y = af["section_offset_y"].as<double>();
            const auto aerodynamic_center = af["aerodynamic_center"].as<double>();
            auto aoa = af["polars"][0]["re_sets"][0]["cl"]["grid"].as<std::vector<double>>();
            std::ranges::transform(aoa, std::begin(aoa), [](auto degrees) {
                return degrees * std::numbers::pi / 180.;
            });
            const auto cl = af["polars"][0]["re_sets"][0]["cl"]["values"].as<std::vector<double>>();
            const auto cd = af["polars"][0]["re_sets"][0]["cd"]["values"].as<std::vector<double>>();
            const auto cm = af["polars"][0]["re_sets"][0]["cm"]["values"].as<std::vector<double>>();

            aero_sections.emplace_back(
                id, s, chord, section_offset_x, section_offset_y, aerodynamic_center, twist, aoa, cl,
                cd, cm
            );
            ++id;
        }

        aero_builder.SetAirfoilSections(0UL, aero_sections);
    }

    //--------------------------------------------------------------------------
    // Interface
    //--------------------------------------------------------------------------

    // Build turbine interface
    auto interface = builder.Build();

    //--------------------------------------------------------------------------
    // Simulation
    //--------------------------------------------------------------------------

    // Loop through solution iterations
    for (auto i : std::views::iota(1U, n_steps)) {
        // Calculate time
        const auto t{static_cast<double>(i) * time_step};

        // Update aero points
        interface.Aerodynamics().CalculateMotion(interface->GetHostState());

        // Set inflow velocity at aero
        for (int nb = 0; nb < interface.Aerodynamics().bodies.size(); ++nb) {
            for (int ns = 0; ns < interface.Aerodynamics().bodies[ns].v_inflow.size(); ++ns) {
                for (int dir = 0; dir < 3; ++dir) {
                    interface.Aerodynamics()
                        .bodies[nb]
                        .v_inflow[ns][dir] = v_inflow[dir];
                }
            }
        }

        // Calculate aero loads
        interface.Aerodynamics().CalculateAerodynamicLoads(fluid_density);

        // Calculate nodal loads
        interface.Aerodynamics().CalculateNodalLoads();

        // Take step
        const auto converged = interface.Step();

        // Check convergence
        ASSERT_EQ(converged, true);
    }
}

}  // namespace kynema::tests
