Example: Rigid body with three springs
======================================

This example will walkthrough how to run simulation of a forced rigid body supported by three nonlinear springs (see :ref:`sec-rb-springs` for details) using Kynema-FMB's high level API.  This model is associated with the DeepCWind platform.
For the most up to date working version of this code, look in ``tests/documentation_tests/floating_platform``.

As with any C++ program, start with the includes.
As a Kokkos-based library, you'll need to include ``Kokkos_Core.hpp`` for setup and teardown operations.
In this example, we will be interfacing to Kynema-FMB using its CFD interface.
In addition to the interface itself, we will include the builder, which is a factory class which will aid in setting up all the necessary data structures.

.. code-block:: cpp

    #include <array>
    #include <cassert>
    #include <Kokkos_Core.hpp>
    #include <interfaces/cfd/interface.hpp>
    #include <interfaces/cfd/interface_builder.hpp>

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

    constexpr auto time_step = 0.01;
    constexpr auto t_end = 1.;      
    constexpr auto rho_inf = 0.0;   
    constexpr auto max_iter = 5;    
    const auto n_steps{static_cast<size_t>(ceil(t_end / time_step)) + 1};
    constexpr auto gravity = std::array{0., 0., -9.8124};
    constexpr auto platform_mass{1.419625E+7}; 
    constexpr auto platform_moi = std::array{1.2898E+10, 1.2851E+10, 1.4189E+10};
    constexpr auto platform_mass_matrix = std::array{
        std::array{platform_mass, 0., 0., 0., 0., 0.},  
        std::array{0., platform_mass, 0., 0., 0., 0.},  
        std::array{0., 0., platform_mass, 0., 0., 0.},  
        std::array{0., 0., 0., platform_moi[0], 0., 0.},
        std::array{0., 0., 0., 0., platform_moi[1], 0.},
        std::array{0., 0., 0., 0., 0., platform_moi[2]},
    };
    constexpr auto mooring_line_stiffness{48.9e3};     
    constexpr auto mooring_line_initial_length{55.432};

We start by creating an interface builder.
You can string together setter functions, or call them individually, depending on what works best with your application.

First, we'll set some general parameters to be used by the solver.
These include the gravity, time step size, numerical damping factor, and maximum number of nonlinear iterations taken per step.

.. code-block:: cpp

    auto interface_builder = kynema_fmb::interfaces::cfd::InterfaceBuilder{}
                              .SetGravity(gravity)
                              .SetTimeStep(time_step)
                              .SetDampingFactor(rho_inf)
                              .SetMaximumNonlinearIterations(max_iter);

The floating platform itself is modeled as a point mass.
To initialize it, we need two pieces of information: its initial position/orientation (as a quaternion) and its mass matrix.

.. code-block:: cpp

    interface_builder.EnableFloatingPlatform(true)
        .SetFloatingPlatformPosition({0., 0., -7.53, 1., 0., 0., 0.})
        .SetFloatingPlatformMassMatrix(platform_mass_matrix);

Mooring lines are modeled as linear springs and require four pieces of information: the stiffness, the initial length, the position of the fairlead point, and the position of the anchor point.
The first argument to each of these setters is the index to the mooring line to be specified.

While the number of Mooring lines must be set before any other information, the other parameters may be set in any order.
We specify all the information for each mooring line at once here, you can also set all the lengths for each line before moving on to the stiffnesses, and so on.
Pick whichever style best fits your application.

.. code-block:: cpp

    interface_builder.SetNumberOfMooringLines(3)
        .SetMooringLineStiffness(0, mooring_line_stiffness)
        .SetMooringLineUndeformedLength(0, mooring_line_initial_length)
        .SetMooringLineFairleadPosition(0, {-40.87, 0.0, -14.})
        .SetMooringLineAnchorPosition(0, {-105.47, 0.0, -58.4})
        .SetMooringLineStiffness(1, mooring_line_stiffness)
        .SetMooringLineUndeformedLength(1, mooring_line_initial_length)
        .SetMooringLineFairleadPosition(1, {20.43, -35.39, -14.})
        .SetMooringLineAnchorPosition(1, {52.73, -91.34, -58.4})
        .SetMooringLineStiffness(2, mooring_line_stiffness)
        .SetMooringLineUndeformedLength(2, mooring_line_initial_length)
        .SetMooringLineFairleadPosition(2, {20.43, 35.39, -14.})
        .SetMooringLineAnchorPosition(2, {52.73, 91.34, -58.4});

When done, call the ``.Build()`` function to generate all of Kynema-FMB's data structures and create the interface itself.

.. code-block:: cpp

   auto interface = interface_builder.Build();

We now compute the buoyancy forces for use during the time stepping process.

..code-block:: cpp

   const auto initial_spring_force = 1907514.4912628897;
   const auto platform_gravity_force = -gravity[2] * platform_mass;
   const auto buoyancy_force = initial_spring_force + platform_gravity_force;

The process of taking each time step is controlled by the user.
Control commands and loads can be changed freely throughout the simulation, either as part of a coupling to an external code or as response to discrete events.

.. code-block:: cpp

   for (auto i = 0U; i < n_steps; ++i) {
   ...
   }

Within this loop, we first set the time-dependent buoyancy forces and moments to the floating platform.

.. code-block:: cpp

   const auto t = static_cast<double>(i) * time_step;
   interface.turbine.floating_platform.node.loads[1] = 1e6 * sin(2. * M_PI / 20. * t);
   interface.turbine.floating_platform.node.loads[2] = buoyancy_force + 0.5 * initial_spring_force * sin(2. * M_PI / 20. * t);
   interface.turbine.floating_platform.node.loads[3] = 5.0e5 * sin(2. * M_PI / 15. * t);
   interface.turbine.floating_platform.node.loads[4] = 1.0e6 * sin(2. * M_PI / 30. * t);
   interface.turbine.floating_platform.node.loads[5] = 2.0e7 * sin(2. * M_PI / 60. * t);

Finally, we call the ``Step`` function to advance the simulation forward one time step.
This function returns a boolean stating if the time step converged or not.
                        
.. code-block:: cpp

   auto converged = interface.Step();
   [[maybe_unused]] const auto converged = interface.Step();
   assert(converged);

And that's it - the simulation will advance the solution in time.
At any time, you can access the current position and orientation of the platform as a quaternion through the ``interface/turbine.floating_platform.node.displacement`` variable.
You can also tell Kynema-FMB to write out the solution to a file at each time step by providing an output file name to the ``InterfaceBuilder::SetOutputFile`` method before building the interface.
