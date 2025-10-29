#include <array>
#include <cassert>

#include <Kokkos_Core.hpp>
#include <elements/beams/hollow_circle_properties.hpp>
#include <interfaces/turbine/turbine_interface.hpp>
#include <interfaces/turbine/turbine_interface_builder.hpp>
#include <yaml-cpp/yaml.h>

int main() {
    // Kynema is based on Kokkos for performance portability.  Make sure to
    // call Kokkos::initialize before creating any Kynema data structures
    // and Kokkos::finalize after all of those data structures have been destroyed.
    Kokkos::initialize();
    {
        const auto duration{0.1};      // Simulation duration in seconds
        const auto time_step{0.01};    // Time step for the simulation
        const auto n_blades{3};        // Number of blades in turbine
        const auto n_blade_nodes{11};  // Number of nodes per blade
        const auto n_tower_nodes{11};  // Number of nodes in tower

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
            .SetMaximumNonlinearIterations(6)
            .SetAbsoluteErrorTolerance(1e-6)
            .SetRelativeErrorTolerance(1e-4);

        // Set output file path
        builder.Outputs().SetOutputFilePath("TurbineInterfaceTest.IEA15");

        // Read WindIO yaml file
        // System parameters will usually be provided in the form of a .yaml file.  Here,
        // we use yaml-cpp to read and parse this file.  Do note that the builder itself does
        // not take in a YAML::Node itself, but requires the user to ferry the information
        // as desired.  Feel free to write your own .yaml parser, use any other file format,
        // or otherwise generate your problem set up as is convienent for your application.
        const auto wio = YAML::LoadFile("./IEA-15-240-RWT.yaml");

        // WindIO components
        // For this case, we will create four parts:
        // 1.  The main tower
        // 2.  The turbine nacelle
        // 3.  The three blades
        // 4.  The turbine hub
        const auto& wio_tower = wio["components"]["tower"];
        const auto& wio_nacelle = wio["components"]["nacelle"];
        const auto& wio_blade = wio["components"]["blade"];
        const auto& wio_hub = wio["components"]["hub"];

        //--------------------------------------------------------------------------
        // Build Turbine
        //--------------------------------------------------------------------------

        // Get turbine builder
        // The Turbine Interface Builder we created above can create a Turbine Builder, which
        // will be used to create the four parts mentioned above, as well as setting the
        // orientation of the various components.
        auto& turbine_builder = builder.Turbine();
        turbine_builder.SetAzimuthAngle(0.)
            .SetRotorApexToHub(0.)
            .SetHubDiameter(wio_hub["diameter"].as<double>())
            .SetConeAngle(wio_hub["cone_angle"].as<double>())
            //.SetBladePitchAngle(M_PI / 2.)
            .SetShaftTiltAngle(wio_nacelle["drivetrain"]["uptilt"].as<double>())
            //.SetNacelleYawAngle(M_PI)
            .SetTowerAxisToRotorApex(wio_nacelle["drivetrain"]["overhang"].as<double>())
            .SetTowerTopToRotorApex(wio_nacelle["drivetrain"]["distance_tt_hub"].as<double>());

        //--------------------------------------------------------------------------
        // Build Blades
        //--------------------------------------------------------------------------

        // Loop through blades and set parameters
        // Each blade is added in reference coordinates and then rotated onto the rotor
        // automatically by Kynema.  The blades of a turbine are assumed to be
        // equally spaced around the rotor - for example, a turbine with three blades
        // will have one blade every 120 degrees about hte X-axis, with the first blade
        // starting along the global Z-axis.
        for (auto j = 0U; j < n_blades; ++j) {
            // Get the blade builder
            // The Turbine Builder will automatically create a new Blade Builder for each
            // blade that is referenced.  The blade numbers should be contiguous, starting
            // at 0, and the blades will be distributed in counter-clockwise order based
            // on their index.  Blades can be added in any order or edited at any time before
            // final assembly by accessing the approprite Blade Builder through this interface.
            auto& blade_builder = turbine_builder.Blade(j);

            // Set blade parameters
            // Kynema's Turbine model assumes one, high order, element per blade.
            // In this case, eleven blade node are used, which we have seen to give accurate
            // results for the IEA15MW blades.
            //
            // The grid points specified in the input will define the number and location of
            // the quadrature points used for assembling the system matrix.
            blade_builder.SetElementOrder(n_blade_nodes - 1).PrescribedRootMotion(false);

            // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
            const auto ref_axis = wio_blade["outer_shape_bem"]["reference_axis"];
            const auto axis_grid = ref_axis["x"]["grid"].as<std::vector<double>>();
            const auto x_values = ref_axis["x"]["values"].as<std::vector<double>>();
            const auto y_values = ref_axis["y"]["values"].as<std::vector<double>>();
            const auto z_values = ref_axis["z"]["values"].as<std::vector<double>>();
            for (auto i = 0U; i < axis_grid.size(); ++i) {
                blade_builder.AddRefAxisPoint(
                    axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                    kynema::interfaces::components::ReferenceAxisOrientation::Z
                );
            }

            // Add reference axis twist
            const auto twist = wio_blade["outer_shape_bem"]["twist"];
            const auto twist_grid = twist["grid"].as<std::vector<double>>();
            const auto twist_values = twist["values"].as<std::vector<double>>();
            for (auto i = 0U; i < twist_grid.size(); ++i) {
                blade_builder.AddRefAxisTwist(twist_grid[i], twist_values[i]);
            }

            // Add blade section properties
            // The stiffness and inertia matrices must be specified for each section.  These are
            // provided to the Blade Builder by way of a std::array<std::array<double, 6>> object.
            // The WindIO file provides them as their tensor components, so we decompress that
            // structure here.
            const auto stiff_matrix =
                wio_blade["elastic_properties_mb"]["six_x_six"]["stiff_matrix"];
            const auto inertia_matrix =
                wio_blade["elastic_properties_mb"]["six_x_six"]["inertia_matrix"];
            const auto k_grid = stiff_matrix["grid"].as<std::vector<double>>();
            const auto m_grid = inertia_matrix["grid"].as<std::vector<double>>();
            const auto n_sections = k_grid.size();
            if (m_grid.size() != k_grid.size()) {
                throw std::runtime_error("stiffness and mass matrices not on same grid");
            }
            for (auto i = 0U; i < n_sections; ++i) {
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
                    kynema::interfaces::components::ReferenceAxisOrientation::Z
                );
            }
        }

        //--------------------------------------------------------------------------
        // Build Tower
        //--------------------------------------------------------------------------

        // Get the tower builder
        // The Turbine Builder can also create a Tower Builder object.  This acts much
        // like the Blade builder class, but there is only one Tower.
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
        for (auto i = 0U; i < axis_grid.size(); ++i) {
            tower_builder.AddRefAxisPoint(
                axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                kynema::interfaces::components::ReferenceAxisOrientation::Z
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
        for (auto i = 0U; i < t_diameter_grid.size(); ++i) {
            // Create section mass and stiffness matrices
            // Kynema provides the helper function GenerateHollowCircleSection to create the
            // inertia and stiffness matrices (M_star and C_star) needed to represent the tower
            const auto section = kynema::beams::GenerateHollowCircleSection(
                t_diameter_grid[i], t_material["E"].as<double>(), t_material["G"].as<double>(),
                t_material["rho"].as<double>(), t_diameter_values[i], t_wall_thickness[i],
                t_material["nu"].as<double>()
            );

            // Add section
            tower_builder.AddSection(
                t_diameter_grid[i], section.M_star, section.C_star,
                kynema::interfaces::components::ReferenceAxisOrientation::Z
            );
        }

        //--------------------------------------------------------------------------
        // Add mass elements
        //--------------------------------------------------------------------------

        // The nacelle and hub are represented as simple mass elements on the system matrix.
        // For these parts, only an inertia matrix is needed.

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

        // Apply load to tower-top node
        interface.Turbine().tower.nodes.back().loads = {1e5, 0., 0., 0., 0., 0.};

        // Apply torque to turbine shaft
        interface.Turbine().torque_control = 1e8;

        // Calculate number of steps
        const auto n_steps{static_cast<size_t>(duration / time_step)};

        // Loop through solution iterations
        // The process of taking each step is controlled by the user.  Control commands
        // and loads can be changed freely throughout the simulation, either as part
        // of a coupling to an external physics code or in respose to discrete events.
        for (auto i = 1U; i < n_steps; ++i) {
            // Calculate time
            const auto t{static_cast<double>(i) * time_step};

            // Set the pitch on blade 3
            interface.Turbine().blade_pitch_control[2] = t * 0.5;

            // Set the yaw angle
            interface.Turbine().yaw_control = t * 0.3;

            // Turn off the torque control after 500 steps
            if (i % 500 == 0) {
                interface.Turbine().torque_control = 0.;
            }

            // Take a single time step
            [[maybe_unused]] const auto converged = interface.Step();

            // Check convergence
            assert(converged);
        }
    }
    // Make sure to call finalize after all Kynema data structures are deleted
    // and you're ready to exit your application.
    Kokkos::finalize();
    return 0;
}
