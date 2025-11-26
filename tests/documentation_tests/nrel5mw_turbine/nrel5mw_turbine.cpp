#include <array>
#include <cassert>

#include <Kokkos_Core.hpp>
#include <elements/beams/hollow_circle_properties.hpp>
#include <interfaces/components/inflow.hpp>
#include <interfaces/turbine/turbine_interface.hpp>
#include <interfaces/turbine/turbine_interface_builder.hpp>
#include <yaml-cpp/yaml.h>

int main() {
    // Kynema is based on Kokkos for performance portability.  Make sure to
    // call Kokkos::initialize before creating any Kynema data structures
    // and Kokkos::finalize after all of those data structures have been destroyed.
    Kokkos::initialize();
    {
        const auto duration{60.};                     // Simulation duration in seconds
        const auto time_step{0.0008608815426997245};  // Time step for the simulation
        const auto n_blades{3};                       // Number of blades in turbine
        const auto n_blade_nodes{11};                 // Number of nodes per blade
        const auto n_tower_nodes{11};                 // Number of nodes in tower
        const auto rotor_speed_init{12.1 * std::numbers::pi / 30.};
        const auto fluid_density{1.0};
        constexpr auto vel_h{11.4};
        constexpr auto h_ref{90.};
        constexpr auto pl_exp{0.};
        constexpr auto flow_angle{30. * std::numbers::pi / 180.};

        // Create interface builder
        // This object is the main interface for building turbines.
        auto builder = kynema::interfaces::TurbineInterfaceBuilder{};

        // Set solution parameters
        // The .Solution() function provides options related to controling the solver,
        // such as the time step, numerical damping factor, and convergence criteria.
        //
        // You can also optionally set an output file in which Kynema will write its
        // solution data each iteration.  If this file is not set, no output will be performed.
        //
        // When using the builder, you can string together setter function calls as seen here,
        // or you can call each of them individually.  The setters are very light weight and
        // will only start constructing the major data structures when .Build() is called on builder.
        builder.Solution()
            .EnableDynamicSolve()
            .SetTimeStep(time_step)
            .SetDampingFactor(0.0)
            .SetGravity({0., 0., -9.81})
            .SetMaximumNonlinearIterations(12)
            .SetAbsoluteErrorTolerance(1e-7)
            .SetRelativeErrorTolerance(1e-6);

        // Set output file path
        builder.Outputs().SetOutputFilePath("SteadyInflow.NREL5MW");

        // Read WindIO yaml file
        // System parameters will usually be provided in the form of a .yaml file.  Here,
        // we use yaml-cpp to read and parse this file.  Do note that the builder itself does
        // not take in a YAML::Node itself, but requires the user to ferry the information
        // as desired.  Feel free to write your own .yaml parser, use any other file format,
        // or otherwise generate your problem set up as is convienent for your application.
        const auto wio = YAML::LoadFile("./NREL-5MW-aero.yaml");

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
            .SetConeAngle(wio_hub["cone_angle"].as<double>() * std::numbers::pi / 180.)
            .SetShaftTiltAngle(
                wio_drivetrain["outer_shape"]["uptilt"].as<double>() * std::numbers::pi / 180.
            )
            .SetTowerAxisToRotorApex(wio_drivetrain["outer_shape"]["overhang"].as<double>())
            .SetTowerTopToRotorApex(wio_drivetrain["outer_shape"]["distance_tt_hub"].as<double>())
            .SetGearBoxRatio(wio_drivetrain["gearbox"]["gear_ratio"].as<double>())
            .SetRotorSpeed(rotor_speed_init)
            .SetNacelleYawAngle(flow_angle);

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
                    kynema::interfaces::components::ReferenceAxisOrientation::Z
                );
            }

            // Add reference axis twist
            const auto blade_twist = wio_blade["outer_shape"]["twist"];
            const auto twist_grid = blade_twist["grid"].as<std::vector<double>>();
            const auto twist_values = blade_twist["values"].as<std::vector<double>>();
            for (auto i : std::views::iota(0U, twist_grid.size())) {
                blade_builder.AddRefAxisTwist(
                    twist_grid[i], -twist_values[i] * std::numbers::pi / 180.
                );
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
                throw std::runtime_error("stiffness and mass matrices not on same grid");
            }
            for (auto i : std::views::iota(0U, n_sections)) {
                if (abs(m_grid[i] - k_grid[i]) > 1e-8) {
                    throw std::runtime_error("stiffness and mass matrices not on same grid");
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
                    kynema::interfaces::components::ReferenceAxisOrientation::Z
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
            .SetElementOrder(n_tower_nodes - 1)  // Set element order to num nodes - 1
            .SetQuadratureStyle(kynema::interfaces::components::BeamInput::QuadratureStyle::Segmented
            )
            .SetSectionRefinement(4)
            .PrescribedRootMotion(false);  // Fix displacement of tower base node

        // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
        const auto t_ref_axis = wio_tower["reference_axis"];
        const auto axis_grid = t_ref_axis["x"]["grid"].as<std::vector<double>>();
        const auto x_values = t_ref_axis["x"]["values"].as<std::vector<double>>();
        const auto y_values = t_ref_axis["y"]["values"].as<std::vector<double>>();
        const auto z_values = t_ref_axis["z"]["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, axis_grid.size())) {
            tower_builder.AddRefAxisPoint(
                axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                kynema::interfaces::components::ReferenceAxisOrientation::Z
            );
        }

        // Set tower base position from first reference axis point
        turbine_builder.SetTowerBasePosition({x_values[0], y_values[0], z_values[0], 1., 0., 0., 0.}
        );

        // Add reference axis twist (zero for tower)
        tower_builder.AddRefAxisTwist(0.0, 0.0).AddRefAxisTwist(1.0, 0.0);

        // Find the tower material properties
        const auto tower_wall_thickness = wio_tower["structure"]["layers"][0]["thickness"];
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
                "Material '" + tower_material_name + "' not found in materials section"
            );
        }

        // Add tower section properties
        const auto elastic_modulus = tower_material["E"].as<double>();
        const auto shear_modulus = tower_material["G"].as<double>();
        const auto poisson_ratio = tower_material["nu"].as<double>();
        const auto density = tower_material["rho"].as<double>();
        const auto tower_diameter = wio_tower["outer_shape"]["outer_diameter"];
        const auto tower_diameter_grid = tower_diameter["grid"].as<std::vector<double>>();
        const auto tower_diameter_values = tower_diameter["values"].as<std::vector<double>>();
        const auto tower_wall_thickness_values =
            tower_wall_thickness["values"].as<std::vector<double>>();

        for (auto i : std::views::iota(0U, tower_diameter_grid.size())) {
            // Create section mass and stiffness matrices
            const auto section = kynema::beams::GenerateHollowCircleSection(
                tower_diameter_grid[i], elastic_modulus, shear_modulus, density,
                tower_diameter_values[i], tower_wall_thickness_values[i], poisson_ratio
            );

            // Add section
            tower_builder.AddSection(
                section.position, section.M_star, section.C_star,
                kynema::interfaces::components::ReferenceAxisOrientation::Z
            );
        }

        //--------------------------------------------------------------------------
        // Add mass elements
        //--------------------------------------------------------------------------

        // Get nacelle mass properties from WindIO
        const auto nacelle_props = wio_drivetrain["elastic_properties"];

        // Set the nacelle inertia matrix in the turbine builder
        if (nacelle_props) {
            // Get the nacelle mass, inertia properties, and location relative to
            // the tower top
            const auto nacelle_mass = nacelle_props["mass"].as<double>();
            const auto nacelle_inertia = nacelle_props["inertia"].as<std::vector<double>>();
            const auto nacelle_cm_offset = nacelle_props["location"].as<std::vector<double>>();

            turbine_builder.SetNacelleInertiaMatrix(
                {{{nacelle_mass, 0., 0., 0., 0., 0.},
                  {0., nacelle_mass, 0., 0., 0., 0.},
                  {0., 0., nacelle_mass, 0., 0., 0.},
                  {0., 0., 0., nacelle_inertia[0], nacelle_inertia[3], nacelle_inertia[4]},
                  {0., 0., 0., nacelle_inertia[3], nacelle_inertia[1], nacelle_inertia[5]},
                  {0., 0., 0., nacelle_inertia[4], nacelle_inertia[5], nacelle_inertia[2]}}},
                {nacelle_cm_offset[0], nacelle_cm_offset[1], nacelle_cm_offset[2]}
            );
        }

        // Set yaw bearing inertia if yaw component exists
        if (wio_yaw) {
            // Get yaw bearing mass properties from WindIO
            const auto yaw_bearing_mass = wio_yaw["elastic_properties"]["mass"].as<double>();
            const auto yaw_bearing_inertia =
                wio_yaw["elastic_properties"]["inertia"].as<std::vector<double>>();

            // Set the yaw bearing inertia matrix in the turbine builder
            turbine_builder.SetYawBearingInertiaMatrix(
                {{{yaw_bearing_mass, 0., 0., 0., 0., 0.},
                  {0., yaw_bearing_mass, 0., 0., 0., 0.},
                  {0., 0., yaw_bearing_mass, 0., 0., 0.},
                  {0., 0., 0., yaw_bearing_inertia[0], 0., 0.},
                  {0., 0., 0., 0., yaw_bearing_inertia[1], 0.},
                  {0., 0., 0., 0., 0., yaw_bearing_inertia[2]}}}
            );
        }

        // Get generator rotational inertia and gearbox ratio from WindIO
        const auto generator_inertia =
            wio_drivetrain["generator"]["elastic_properties"]["inertia"].as<std::vector<double>>();
        const auto gearbox_ratio = wio_drivetrain["gearbox"]["gear_ratio"].as<double>();

        // Get hub mass properties from WindIO
        const auto hub_mass = wio_hub["elastic_properties"]["mass"].as<double>();
        const auto hub_inertia = wio_hub["elastic_properties"]["inertia"].as<std::vector<double>>();

        // Set the hub inertia matrix in the turbine builder
        turbine_builder.SetHubInertiaMatrix(
            {{{hub_mass, 0., 0., 0., 0., 0.},
              {0., hub_mass, 0., 0., 0., 0.},
              {0., 0., hub_mass, 0., 0., 0.},
              {0., 0., 0., hub_inertia[0] + generator_inertia[0] * gearbox_ratio * gearbox_ratio,
               hub_inertia[3], hub_inertia[4]},
              {0., 0., 0., hub_inertia[3], hub_inertia[1], hub_inertia[5]},
              {0., 0., 0., hub_inertia[4], hub_inertia[5], hub_inertia[2]}}}
        );

        //--------------------------------------------------------------------------
        // Build Aerodynamics
        //--------------------------------------------------------------------------

        auto& aero_builder =
            builder.Aerodynamics().EnableAero().SetNumberOfAirfoils(1UL).SetAirfoilToBladeMap(
                std::array{0UL, 0UL, 0UL}
            );

        const auto& airfoil_io = wio["airfoils"];
        auto aero_sections = std::vector<kynema::interfaces::components::AerodynamicSection>{};
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

        //--------------------------------------------------------------------------
        // Build Controller
        //--------------------------------------------------------------------------

        const auto controller_shared_lib_path =
            std::string{static_cast<const char*>(Kynema_ROSCO_LIBRARY)};
        const auto controller_function_name = std::string{"DISCON"};
        const auto controller_input_file = std::string{"DISCON.IN"};
        const auto controller_output_file = std::string{"controller_NREL5MW"};

        auto controller_builder = builder.Controller()
                                      .SetLibraryPath(controller_shared_lib_path)
                                      .SetFunctionName(controller_function_name)
                                      .SetInputFilePath(controller_input_file)
                                      .SetControllerInput(controller_output_file);

        //--------------------------------------------------------------------------
        // Interface
        //--------------------------------------------------------------------------

        // Build turbine interface
        // Now that we are done setting up the system, call build on the initial Inerface
        // Builder that we made back at the beginning.  This step is where everything
        // actually gets sized and allocated for solving the system.  You interact with
        // this system through the Turbine Interface object that's created.
        auto interface = builder.Build();

        //--------------------------------------------------------------------------
        // Simulation
        //--------------------------------------------------------------------------

        // Calculate number of steps
        const auto n_steps{static_cast<size_t>(duration / time_step)};

        auto inflow =
            kynema::interfaces::components::Inflow::SteadyWind(vel_h, h_ref, pl_exp, flow_angle);

        const bool might_crash{false};
        if (might_crash) {
            interface.CloseOutputFile();
        }

        // Loop through solution iterations
        // The process of taking each step is controlled by the user.  Control commands
        // and loads can be changed freely throughout the simulation, either as part
        // of a coupling to an external physics code or in respose to discrete events.
        for (auto i = 1U; i <= n_steps; ++i) {
            // Calculate time
            const auto t{static_cast<double>(i) * time_step};

            std::cout << "Step " << i << ", t = " << std::fixed << std::setprecision(6) << t << " s "
                      << " | " << std::setprecision(2) << t / duration * 100. << "%\n";

            interface.UpdateAerodynamicLoads(fluid_density, [t, &inflow](const std::array<double, 3>& pos) {
                    return inflow.Velocity(t, pos);
                }
            );

            interface.ApplyController(t);

            interface.Turbine().yaw_control = flow_angle;

            // Take a single time step
            [[maybe_unused]] const auto converged = interface.Step();

            // Check convergence
            assert(converged);

            // Do output every 20 steps
            if (i % 20 == 0) {
                std::cout << " --- Writing output file\n";
                if (might_crash) {
                    interface.OpenOutputFile();
                }
                interface.WriteOutput();
                if (might_crash) {
                    interface.CloseOutputFile();
                }
            }
        }
    }
    // Make sure to call finalize after all Kynema data structures are deleted
    // and you're ready to exit your application.
    Kokkos::finalize();
    return 0;
}
