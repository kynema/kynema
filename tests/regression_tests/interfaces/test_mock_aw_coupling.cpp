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
    constexpr auto rotor_speed_init{0.};
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
    const auto& wio_blade = wio["components"]["blade"];
    const auto& wio_tower = wio["components"]["tower"];
    const auto& wio_drivetrain = wio["components"]["drivetrain"];
    const auto& wio_hub = wio["components"]["hub"];
    const auto& wio_yaw = wio["components"]["yaw"];

    //--------------------------------------------------------------------------
    // Build Turbine
    //--------------------------------------------------------------------------

    // Get turbine builder
    auto& turbine_builder = builder.Turbine();
    turbine_builder.SetAzimuthAngle(0.)
        .SetRotorApexToHub(0.)
        .SetHubDiameter(wio_hub["diameter"].as<double>())
        .SetConeAngle(
            wio_hub["cone_angle"].as<double>() * std::numbers::pi / 180.)
        .SetShaftTiltAngle(
            wio_drivetrain["outer_shape"]["uptilt"].as<double>() *
            std::numbers::pi / 180.)
        .SetTowerAxisToRotorApex(
            wio_drivetrain["outer_shape"]["overhang"].as<double>())
        .SetTowerTopToRotorApex(
            wio_drivetrain["outer_shape"]["distance_tt_hub"].as<double>())
        .SetGearBoxRatio(wio_drivetrain["gearbox"]["gear_ratio"].as<double>())
        .SetRotorSpeed(rotor_speed_init);

    //--------------------------------------------------------------------------
    // Build Blades
    //--------------------------------------------------------------------------

    // Loop through blades and set parameters
    for (auto j : std::views::iota(0U, static_cast<uint>(n_blades))) {
        // Get the blade builder
        auto& blade_builder = turbine_builder.Blade(j);

        // Set blade parameters
        blade_builder.SetElementOrder(n_blade_nodes - 1)
            .PrescribedRootMotion(false)
            .SetSectionRefinement(2);

        // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
        const auto ref_axis = wio_blade["reference_axis"];
        const auto axis_grid = ref_axis["x"]["grid"].as<std::vector<double>>();
        const auto x_values = ref_axis["x"]["values"].as<std::vector<double>>();
        const auto y_values = ref_axis["y"]["values"].as<std::vector<double>>();
        const auto z_values = ref_axis["z"]["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, axis_grid.size())) {
            blade_builder.AddRefAxisPoint(
                axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                kynema::interfaces::components::ReferenceAxisOrientation::Z);
        }

        // Add reference axis twist
        const auto blade_twist = wio_blade["outer_shape"]["twist"];
        const auto twist_grid = blade_twist["grid"].as<std::vector<double>>();
        const auto twist_values =
            blade_twist["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, twist_grid.size())) {
            blade_builder.AddRefAxisTwist(
                twist_grid[i], -twist_values[i] * std::numbers::pi / 180.);
        }

        const auto inertia_matrix =
            wio_blade["structure"]["elastic_properties"]["inertia_matrix"];
        const auto stiffness_matrix =
            wio_blade["structure"]["elastic_properties"]["stiffness_matrix"];

        // Add blade section properties
        const auto k_grid = stiffness_matrix["grid"].as<std::vector<double>>();
        const auto m_grid = inertia_matrix["grid"].as<std::vector<double>>();
        const auto n_sections = k_grid.size();
        if (m_grid.size() != k_grid.size()) {
            throw std::runtime_error(
                "stiffness and mass matrices not on same grid");
        }
        for (auto i : std::views::iota(0U, n_sections)) {
            if (abs(m_grid[i] - k_grid[i]) > 1e-8) {
                throw std::runtime_error(
                    "stiffness and mass matrices not on same grid");
            }
            const auto mass = inertia_matrix["mass"][i].as<double>();
            const auto cm_x = inertia_matrix["cm_x"][i].as<double>();
            const auto cm_y = inertia_matrix["cm_y"][i].as<double>();
            const auto i_cp = inertia_matrix["i_cp"][i].as<double>();
            const auto i_edge = inertia_matrix["i_edge"][i].as<double>();
            const auto i_flap = inertia_matrix["i_flap"][i].as<double>();
            const auto i_plr = inertia_matrix["i_plr"][i].as<double>();

            const auto k11 = stiffness_matrix["K11"][i].as<double>();
            const auto k12 = stiffness_matrix["K12"][i].as<double>();
            const auto k13 = stiffness_matrix["K13"][i].as<double>();
            const auto k14 = stiffness_matrix["K14"][i].as<double>();
            const auto k15 = stiffness_matrix["K15"][i].as<double>();
            const auto k16 = stiffness_matrix["K16"][i].as<double>();
            const auto k22 = stiffness_matrix["K22"][i].as<double>();
            const auto k23 = stiffness_matrix["K23"][i].as<double>();
            const auto k24 = stiffness_matrix["K24"][i].as<double>();
            const auto k25 = stiffness_matrix["K25"][i].as<double>();
            const auto k26 = stiffness_matrix["K26"][i].as<double>();
            const auto k33 = stiffness_matrix["K33"][i].as<double>();
            const auto k34 = stiffness_matrix["K34"][i].as<double>();
            const auto k35 = stiffness_matrix["K35"][i].as<double>();
            const auto k36 = stiffness_matrix["K36"][i].as<double>();
            const auto k44 = stiffness_matrix["K44"][i].as<double>();
            const auto k45 = stiffness_matrix["K45"][i].as<double>();
            const auto k46 = stiffness_matrix["K46"][i].as<double>();
            const auto k55 = stiffness_matrix["K55"][i].as<double>();
            const auto k56 = stiffness_matrix["K56"][i].as<double>();
            const auto k66 = stiffness_matrix["K66"][i].as<double>();

            blade_builder.AddSection(
                m_grid[i],
                {{
                    {mass, 0., 0., 0., 0., -mass * cm_y},
                    {0., mass, 0., 0., 0., mass * cm_x},
                    {0., 0., mass, mass * cm_y, -mass * cm_x, 0.},
                    {0., 0., mass * cm_y, i_edge, -i_cp, 0.},
                    {0., 0., -mass * cm_x, -i_cp, i_flap, 0.},
                    {-mass * cm_y, mass * cm_x, 0., 0., 0., i_plr},
                }},
                {{
                    {k11, k12, k13, k14, k15, k16},
                    {k12, k22, k23, k24, k25, k26},
                    {k13, k23, k33, k34, k35, k36},
                    {k14, k24, k34, k44, k45, k46},
                    {k15, k25, k35, k45, k55, k56},
                    {k16, k26, k36, k46, k56, k66},
                }},
                kynema::interfaces::components::ReferenceAxisOrientation::Z);
        }
    }

    //--------------------------------------------------------------------------
    // Build Tower
    //--------------------------------------------------------------------------

    // Get the tower builder
    auto& tower_builder = turbine_builder.Tower();

    // Set tower parameters
    tower_builder
        .SetElementOrder(
            n_tower_nodes - 1)        // Set element order to num nodes - 1
        .PrescribedRootMotion(false); // Fix displacement of tower base node

    // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
    const auto t_ref_axis = wio_tower["reference_axis"];
    const auto axis_grid = t_ref_axis["x"]["grid"].as<std::vector<double>>();
    const auto x_values = t_ref_axis["x"]["values"].as<std::vector<double>>();
    const auto y_values = t_ref_axis["y"]["values"].as<std::vector<double>>();
    const auto z_values = t_ref_axis["z"]["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, axis_grid.size())) {
        tower_builder.AddRefAxisPoint(
            axis_grid[i], {x_values[i], y_values[i], z_values[i]},
            kynema::interfaces::components::ReferenceAxisOrientation::Z);
    }

    // Set tower base position from first reference axis point
    turbine_builder.SetTowerBasePosition(
        {x_values[0], y_values[0], z_values[0], 1., 0., 0., 0.});

    // Add reference axis twist (zero for tower)
    tower_builder.AddRefAxisTwist(0.0, 0.0).AddRefAxisTwist(1.0, 0.0);

    // Find the tower material properties
    const auto tower_diameter = wio_tower["outer_shape"]["outer_diameter"];
    const auto tower_wall_thickness =
        wio_tower["structure"]["layers"][0]["thickness"];
    const auto tower_material_name =
        wio_tower["structure"]["layers"][0]["material"].as<std::string>();

    YAML::Node tower_material;
    for (const auto& m : wio["materials"]) {
        if (m["name"] && m["name"].as<std::string>() == tower_material_name) {
            tower_material = m.as<YAML::Node>();
            break;
        }
    }
    if (!tower_material) {
        throw std::runtime_error(
            "Material '" + tower_material_name +
            "' not found in materials section");
    }

    // Add tower section properties
    const auto elastic_modulus = tower_material["E"].as<double>();
    const auto shear_modulus = tower_material["G"].as<double>();
    const auto poisson_ratio = tower_material["nu"].as<double>();
    const auto density = tower_material["rho"].as<double>();
    for (auto i : std::views::iota(0U, tower_diameter["grid"].size())) {
        // Create section mass and stiffness matrices
        const auto section = kynema::beams::GenerateHollowCircleSection(
            tower_diameter["grid"][i].as<double>(), elastic_modulus,
            shear_modulus, density, tower_diameter["values"][i].as<double>(),
            tower_wall_thickness["values"][i].as<double>(), poisson_ratio);

        // Add section
        tower_builder.AddSection(
            tower_diameter["grid"][i].as<double>(), section.M_star,
            section.C_star,
            kynema::interfaces::components::ReferenceAxisOrientation::Z);
    }

    //--------------------------------------------------------------------------
    // Add mass elements
    //--------------------------------------------------------------------------

    // Get nacelle mass properties from WindIO
    const auto nacelle_props = wio_drivetrain["elastic_properties"];
    const auto nacelle_mass = nacelle_props["mass"].as<double>();
    const auto nacelle_inertia =
        nacelle_props["inertia"].as<std::vector<double>>();

    // Nacelle center of mass offset from yaw bearing
    const auto nacelle_cm_offset =
        nacelle_props["location"].as<std::vector<double>>();

    // Set the nacelle inertia matrix in the turbine builder
    turbine_builder.SetNacelleInertiaMatrix(
        {{{nacelle_mass, 0., 0., 0., 0., 0.},
          {0., nacelle_mass, 0., 0., 0., 0.},
          {0., 0., nacelle_mass, 0., 0., 0.},
          {0., 0., 0., nacelle_inertia[0], nacelle_inertia[3],
           nacelle_inertia[4]},
          {0., 0., 0., nacelle_inertia[3], nacelle_inertia[1],
           nacelle_inertia[5]},
          {0., 0., 0., nacelle_inertia[4], nacelle_inertia[5],
           nacelle_inertia[2]}}},
        {nacelle_cm_offset[0], nacelle_cm_offset[1], nacelle_cm_offset[2]});

    // Get yaw bearing mass properties from WindIO
    const auto yaw_bearing_mass =
        wio_yaw["elastic_properties"]["mass"].as<double>();
    const auto yaw_bearing_inertia =
        wio_yaw["elastic_properties"]["inertia"].as<std::vector<double>>();

    // Set the yaw bearing inertia matrix in the turbine builder
    turbine_builder.SetYawBearingInertiaMatrix(
        {{{yaw_bearing_mass, 0., 0., 0., 0., 0.},
          {0., yaw_bearing_mass, 0., 0., 0., 0.},
          {0., 0., yaw_bearing_mass, 0., 0., 0.},
          {0., 0., 0., yaw_bearing_inertia[0], 0., 0.},
          {0., 0., 0., 0., yaw_bearing_inertia[1], 0.},
          {0., 0., 0., 0., 0., yaw_bearing_inertia[2]}}});

    // Get generator rotational inertia and gearbox ratio from WindIO
    const auto generator_inertia =
        wio_drivetrain["generator"]["elastic_properties"]["inertia"]
            .as<std::vector<double>>();
    const auto gearbox_ratio =
        wio_drivetrain["gearbox"]["gear_ratio"].as<double>();

    // Get hub mass properties from WindIO
    const auto hub_mass = wio_hub["elastic_properties"]["mass"].as<double>();
    const auto hub_inertia =
        wio_hub["elastic_properties"]["inertia"].as<std::vector<double>>();

    // Set the hub inertia matrix in the turbine builder
    turbine_builder.SetHubInertiaMatrix(
        {{{hub_mass, 0., 0., 0., 0., 0.},
          {0., hub_mass, 0., 0., 0., 0.},
          {0., 0., hub_mass, 0., 0., 0.},
          {0., 0., 0., hub_inertia[0] + generator_inertia[0] * gearbox_ratio,
           hub_inertia[3], hub_inertia[4]},
          {0., 0., 0., hub_inertia[3], hub_inertia[1], hub_inertia[5]},
          {0., 0., 0., hub_inertia[4], hub_inertia[5], hub_inertia[2]}}});

    //--------------------------------------------------------------------------
    // Build Aerodynamics
    //--------------------------------------------------------------------------

    auto& aero_builder = builder.Aerodynamics()
                             .EnableAero()
                             .SetNumberOfAirfoils(1UL)
                             .SetAirfoilToBladeMap(std::array{0UL, 0UL, 0UL});

    const auto& airfoil_io = wio["airfoils"];
    auto aero_sections =
        std::vector<kynema::interfaces::components::AerodynamicSection>{};
    auto id = 0UL;
    for (const auto& af : airfoil_io) {
        const auto s = af["spanwise_position"].as<double>();
        const auto chord = af["chord"].as<double>();
        const auto twist = af["twist"].as<double>() * std::numbers::pi / 180.;
        const auto section_offset_x = af["section_offset_x"].as<double>();
        const auto section_offset_y = af["section_offset_y"].as<double>();
        const auto aerodynamic_center = af["aerodynamic_center"].as<double>();
        auto aoa = af["polars"][0]["re_sets"][0]["cl"]["grid"]
                       .as<std::vector<double>>();
        std::ranges::transform(aoa, std::begin(aoa), [](auto degrees) {
            return degrees * std::numbers::pi / 180.;
        });
        const auto cl = af["polars"][0]["re_sets"][0]["cl"]["values"]
                            .as<std::vector<double>>();
        const auto cd = af["polars"][0]["re_sets"][0]["cd"]["values"]
                            .as<std::vector<double>>();
        const auto cm = af["polars"][0]["re_sets"][0]["cm"]["values"]
                            .as<std::vector<double>>();
        aero_sections.emplace_back(
            id, s, chord, section_offset_x, section_offset_y,
            aerodynamic_center, twist, aoa, cl, cd, cm);
        ++id;
    }

    aero_builder.SetAirfoilSections(0UL, aero_sections);

    //--------------------------------------------------------------------------
    // Interface
    //--------------------------------------------------------------------------

    // Build turbine interface
    auto interface = builder.Build();

    //--------------------------------------------------------------------------
    // Simulation
    //--------------------------------------------------------------------------

    // Loop through solution iterations
    for (auto i : std::views::iota(0U, n_steps)) {
        // Calculate time
        const auto t{static_cast<double>(i) * time_step};

        // Update aero points
        interface.Aerodynamics().CalculateMotion(interface.GetHostState());

        // Set inflow velocity at aero
        for (uint nb = 0; nb < interface.Aerodynamics().bodies.size(); ++nb) {
            for (uint ns = 0; ns < interface.Aerodynamics().bodies[nb].v_inflow.size(); ++ns) {
                for (uint dir = 0; dir < 3U; ++dir) {
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
