Example: IEA15MW Turbine
========================

This example will walkthrough how to simulate a wind turbine using Kynema's high level API.
For the most up to date working version of this code, look in ``tests/documentation_tests/iea15mw_turbine``.

As with any C++ program, start with the includes.
We will need to include ``Kokkos_Core.hpp`` to initialize and finalize Kokkos at the beginning of the application. 
We will also include the turbine interface, which we will be using in this example, and its builder.
We'll also include the ``HollowCircleProperties.hpp`` header to support creating the hub of our turbine.
Finally, we'll include the ``yaml-cpp`` file for reading a Wind-IO input file.

.. code-block:: cpp

    #include <array>
    #include <cassert>
    #include <Kokkos_Core.hpp>
    #include <elements/beams/hollow_circle_properties.hpp>
    #include <interfaces/turbine/turbine_interface.hpp>
    #include <interfaces/turbine/turbine_interface_builder.hpp>
    #include <yaml-cpp/yaml.h>

We now create the main function of our program and initialize Kokkos.
We also start a scope around the rest of the program.
This scoping is necessary to ensure that all Kokkos objects are destroyed before finalize is called and the program exits.
Failure to do this will result in lots of nasty, hard to decipher errors on program termination.

.. code-block:: cpp

    int main() {
        Kokkos::initialize();
        {
            ...
        }
        Kokkos::finalize();
        return 0;
    }

Now set up some constants that we'll need when running the simulation.
In general, these will be read fom an input file, but we just define them inline here.

.. code-block:: cpp

    const auto duration{0.1};    
    const auto time_step{0.01};  
    const auto n_blades{3};      
    const auto n_blade_nodes{11};
    const auto n_tower_nodes{11};

Next, we create the turbine interface builder.
This is a factory class which will aid us in setting up all of the structures Kynema needs for its simulation.

.. code-block:: cpp

    auto builder = kynema::interfaces::TurbineInterfaceBuilder{};

The builder has several groupings of options that we will use to set up our problem.
The ``.Solution()`` function provides options related to controling the solver, such as the time step, numerical damping, and convergence criteria.
You can also optionally set an output file in which Kynema will write its solution data each iteration.
If this file is not set, no output will be performed.

When using the builder, you can string together setter function calls as seen here, or you can call each of them individually.
The setters are very light weight and will only start constructing the major data structures when .Build() is called on builder.

.. code-block:: cpp

   builder.Solution()
       .EnableDynamicSolve()
       .SetTimeStep(time_step)
       .SetDampingFactor(0.0)
       .SetGravity({0., 0., -9.81})
       .SetMaximumNonlinearIterations(6)
       .SetAbsoluteErrorTolerance(1e-6)
       .SetRelativeErrorTolerance(1e-4)
       .SetOutputFile("TurbineInterfaceTest.IEA15");

System parameters will usually be provided in the form of a .yaml file.
Here, we use yaml-cpp to read and parse this file.
Do note that the builder itself does not take in a YAML::Node itself, but requires the user to ferry the information as desired.
Feel free to write your own .yaml parser, use any other file format, or otherwise generate your problem set up as is convienent for your application.

.. code-block:: cpp

   const auto wio = YAML::LoadFile("./IEA-15-240-RWT.yaml");

For this problem, we will create four parts: the main tower, the turbine nacelle, the three blades, and the turbine hub

.. code-block:: cpp

    const auto& wio_tower = wio["components"]["tower"];
    const auto& wio_nacelle = wio["components"]["nacelle"];
    const auto& wio_blade = wio["components"]["blade"];
    const auto& wio_hub = wio["components"]["hub"];

The Turbine Interface Builder we created above can create a Turbine Builder, which will be used to create the four parts mentioned above, as well as setting the orientation of the various components.

.. code-block:: cpp

    auto& turbine_builder = builder.Turbine();
    turbine_builder.SetAzimuthAngle(0.)
        .SetRotorApexToHub(0.)
        .SetHubDiameter(wio_hub["diameter"].as<double>())
        .SetConeAngle(wio_hub["cone_angle"].as<double>())
        .SetShaftTiltAngle(wio_nacelle["drivetrain"]["uptilt"].as<double>())
        .SetTowerAxisToRotorApex(wio_nacelle["drivetrain"]["overhang"].as<double>())
        .SetTowerTopToRotorApex(wio_nacelle["drivetrain"]["distance_tt_hub"].as<double>());

Each blade is added in reference coordinates and then rotated onto the rotor automatically by Kynema.
The blades of a turbine are assumed to be equally spaced around the rotor - for example, a turbine with three blades will have one blade every 120 degrees about the X-axis, with the first blade starting along the global Z-axis.

.. code-block:: cpp

    for (auto j = 0U; j < n_blades; ++j) {
    ...
    }

The Turbine Builder will automatically create a new Blade Builder for each blade that is referenced.
The blade numbers should be contiguous, starting at 0, and the blades will be distributed in counter-clockwise order based on their index.
Blades can be added in any order or edited at any time before final assembly by accessing the approprite Blade Builder through this interface.

.. code-block:: cpp

   auto& blade_builder = turbine_builder.Blade(j);

Kynema's Turbine model assumes one, high order, element per blade.
In this case, eleven blade node are used, which we have seen to give accurate results for the IEA15MW blades.

The grid points specified in the input will define the number and location of the quadrature points used for assembling the system matrix.

.. code-block:: cpp

    blade_builder.SetElementOrder(n_blade_nodes - 1).PrescribedRootMotion(false);

We now add the reference axis coordinates of the nodes along the blade.
Note that Wind-IO uses the Z-axis as its reference axis, but your application may differ in this choice.

.. code-block:: cpp

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

Next, we add blade twist about the reference axis.

.. code-block:: cpp

    const auto twist = wio_blade["outer_shape_bem"]["twist"];
    const auto twist_grid = twist["grid"].as<std::vector<double>>();
    const auto twist_values = twist["values"].as<std::vector<double>>();
    for (auto i = 0U; i < twist_grid.size(); ++i) {
        blade_builder.AddRefAxisTwist(twist_grid[i], twist_values[i]);
    }

The stiffness and inertia matrices must be specified for each section.
These are provided to the Blade Builder by way of a std::array<std::array<double, 6>> object.
The WindIO file provides them as their tensor components, so we decompress that structure here.

.. code-block:: cpp

    const auto stiff_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["stiff_matrix"];
    const auto inertia_matrix = wio_blade["elastic_properties_mb"]["six_x_six"]["inertia_matrix"];
    const auto k_grid = stiff_matrix["grid"].as<std::vector<double>>();
    const auto m_grid = inertia_matrix["grid"].as<std::vector<double>>();
    const auto n_sections = k_grid.size();
    for (auto i = 0U; i < n_sections; ++i) {
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

With the blades of the turbine built, the Turbine Builder can also create a Tower Builder object.
This acts much like the Blade builder class, but there is only one Tower.

.. code-block:: cpp
   
    auto& tower_builder = turbine_builder.Tower();
    tower_builder.SetElementOrder(n_tower_nodes - 1).PrescribedRootMotion(false);
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
    const auto tower_base_position =
        std::array{x_values[0], y_values[0], z_values[0], 1., 0., 0., 0.};
    turbine_builder.SetTowerBasePosition(tower_base_position);
    tower_builder.AddRefAxisTwist(0.0, 0.0).AddRefAxisTwist(1.0, 0.0);
    const auto t_layer = wio_tower["internal_structure_2d_fem"]["layers"][0];
    const auto t_material_name = t_layer["material"].as<std::string>();
    YAML::Node t_material;
    for (const auto& m : wio["materials"]) {
        if (m["name"] && m["name"].as<std::string>() == t_material_name) {
            t_material = m.as<YAML::Node>();
            break;
        }
    }
    const auto t_diameter = wio_tower["outer_shape_bem"]["outer_diameter"];
    const auto t_diameter_grid = t_diameter["grid"].as<std::vector<double>>();
    const auto t_diameter_values = t_diameter["values"].as<std::vector<double>>();
    const auto t_wall_thickness = t_layer["thickness"]["values"].as<std::vector<double>>();
   for (auto i = 0U; i < t_diameter_grid.size(); ++i) {
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

The nacelle and hub are represented as simple mass elements on the system matrix.
For these parts, only an inertia matrix is needed.

.. code-block:: cpp
   
    const auto& nacelle_props = wio_nacelle["elastic_properties_mb"];
    const auto system_mass = nacelle_props["system_mass"].as<double>();
    const auto yaw_mass = nacelle_props["yaw_mass"].as<double>();
    const auto system_inertia_tt = nacelle_props["system_inertia_tt"].as<std::vector<double>>();
    const auto total_mass = system_mass + yaw_mass;
    const auto nacelle_inertia_matrix = std::array<std::array<double, 6>, 6>{
        {{total_mass, 0., 0., 0., 0., 0.},
         {0., total_mass, 0., 0., 0., 0.},
         {0., 0., total_mass, 0., 0., 0.},
         {0., 0., 0., system_inertia_tt[0], system_inertia_tt[3], system_inertia_tt[4]},
         {0., 0., 0., system_inertia_tt[3], system_inertia_tt[1], system_inertia_tt[5]},
         {0., 0., 0., system_inertia_tt[4], system_inertia_tt[5], system_inertia_tt[2]}}
    };
    turbine_builder.SetYawBearingInertiaMatrix(nacelle_inertia_matrix);
    const auto& hub_props = wio_hub["elastic_properties_mb"];
    const auto hub_mass = hub_props["system_mass"].as<double>();
    const auto hub_inertia = hub_props["system_inertia"].as<std::vector<double>>();
    const auto hub_inertia_matrix = std::array<std::array<double, 6>, 6>{
        {{hub_mass, 0., 0., 0., 0., 0.},
         {0., hub_mass, 0., 0., 0., 0.},
         {0., 0., hub_mass, 0., 0., 0.},
         {0., 0., 0., hub_inertia[0], hub_inertia[3], hub_inertia[4]},
         {0., 0., 0., hub_inertia[3], hub_inertia[1], hub_inertia[5]},
         {0., 0., 0., hub_inertia[4], hub_inertia[5], hub_inertia[2]}}
    };
    turbine_builder.SetHubInertiaMatrix(hub_inertia_matrix);

Now that we are done setting up the system, call build on the initial Inerface Builder that we made back at the beginning.
This step is where everything actually gets sized and allocated for solving the system.
You interact with this system through the Turbine Interface object that's created.

.. code-block:: cpp

    auto interface = builder.Build();

We now set the initial loads and torque on the turbine.

.. code-block:: cpp

    interface.Turbine().tower.nodes.back().loads = {1e5, 0., 0., 0., 0., 0.};
    interface.Turbine().torque_control = 1e8;

The process of taking each time step is controlled by the user.
Control commands and loads can be changed freely throughout the simulation, either as part of a coupling to an external code or as response to discrete events.

.. code-block:: cpp

    const auto n_steps{static_cast<size_t>(duration / time_step)};
    for (auto i = 1U; i < n_steps; ++i) {
    ...
    }

Within each time step, we set the control commands.

.. code-block:: cpp

    const auto t{static_cast<double>(i) * time_step};
    interface.Turbine().blade_pitch_control[2] = t * 0.5;
    interface.Turbine().yaw_control = t * 0.3;
    if (i % 500 == 0) {
        interface.Turbine().torque_control = 0.;
    }

Finally, we call the ``Step()`` function on the turbine interface, which actually performs the action of advancing the solution in time.
This function returns a boolean stating if Kynema's solver converged or not, which can be checked for error handling.

.. code-block:: cpp

    [[maybe_unused]] const auto converged = interface.Step();
    assert(converged);

And that's it, the simulation will advance to the total solution time.
Kynema will have written out the solution at each time step in NetCDF format.
At any time, the solution can be accessed by looking at the ``interface.Turbine()`` object.
For example, the position and orientation quaternion can be accessed by calling ``interface.Turbine().tower.nodes.back()``.

