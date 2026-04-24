#include <array>
#include <cassert>

#include <Kokkos_Core.hpp>
#include <interfaces/cfd/interface.hpp>
#include <interfaces/cfd/interface_builder.hpp>

int main() {
    // Kynema is based on Kokkos for performance portability.  Make sure to
    // call Kokkos::initialize before creating any Kynema data structures
    // and Kokkos::finalize after all of those data structures have been destroyed.
    Kokkos::initialize();
    {
        // Solution parameters
        constexpr auto time_step = 0.01;  // Time step size
        constexpr auto t_end = 1.;        // Final time of simulation
        constexpr auto rho_inf = 0.0;     // Time stepping damping factor
        constexpr auto max_iter = 5;      // Maximum number of nonlinear steps per time ste[
        const auto n_steps{
            static_cast<size_t>(ceil(t_end / time_step)) + 1
        };  // Number of time steps

        // Define gravity vector
        constexpr auto gravity = std::array{0., 0., -9.8124};  // m/s/s

        // Construct platform mass matrix as a 6x6 "array of arrays"
        constexpr auto platform_mass{1.419625E+7};                                     // kg
        constexpr auto platform_moi = std::array{1.2898E+10, 1.2851E+10, 1.4189E+10};  // kg*m*m
        constexpr auto platform_mass_matrix = std::array{
            std::array{platform_mass, 0., 0., 0., 0., 0.},    // Row 1
            std::array{0., platform_mass, 0., 0., 0., 0.},    // Row 2
            std::array{0., 0., platform_mass, 0., 0., 0.},    // Row 3
            std::array{0., 0., 0., platform_moi[0], 0., 0.},  // Row 4
            std::array{0., 0., 0., 0., platform_moi[1], 0.},  // Row 5
            std::array{0., 0., 0., 0., 0., platform_moi[2]},  // Row 6
        };

        // Define mooring line stiffness and initial length
        constexpr auto mooring_line_stiffness{48.9e3};       // N
        constexpr auto mooring_line_initial_length{55.432};  // m

        // Create cfd interface
        // This interface is the primary way for interacting with Kynema
        // and solving floating platform problems.  You can string together setter
        // functions, or call them individually, depending on which
        // works best with your applications.
        //
        // First, general solution parameters should be set - the gravity force,
        // time step, numerical damping factor, and maximum number of nonlinear
        // iterations.
        auto interface_builder = kynema_fmb::interfaces::cfd::InterfaceBuilder{}
                                     .SetGravity(gravity)
                                     .SetTimeStep(time_step)
                                     .SetDampingFactor(rho_inf)
                                     .SetMaximumNonlinearIterations(max_iter);

        // The floating platfrom is modeled as a point mass and needs two pieces of information:
        // its initial position (as a quaternion) and its mass matrix with inertia information.
        interface_builder.EnableFloatingPlatform(true)
            .SetFloatingPlatformPosition({0., 0., -7.53, 1., 0., 0., 0.})
            .SetFloatingPlatformMassMatrix(platform_mass_matrix);

        // Mooring lines are modeled as springs require three pieces of informations: the stiffness,
        // initial length, and position of each the fairlead and anchor point.  The first
        // argument to each of these setters is the index of the mooring line to be specified.
        //
        // While the number of mooring lines in the problem must be set before other information,
        // the other parameters may be set in any order.  While we specify all of the
        // parameters for each mooring line at once here, you may specify each line's stiffness
        // before specifying their initial lengths, then node positions.  Pick whichever style
        // best fits your application.
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

        // Call the .Build() function when done to generate the interface and all of Kynema's
        // data structures.
        auto interface = interface_builder.Build();

        // Calculate buoyancy force as percentage of gravitational force plus spring forces times
        const auto initial_spring_force = 1907514.4912628897;
        const auto platform_gravity_force = -gravity[2] * platform_mass;
        const auto buoyancy_force = initial_spring_force + platform_gravity_force;

        // Iterate through time steps
        // The process of taking each step is controlled by the user.  Control commands
        // and loads can be changed freely throughout the simulation, either as part
        // of a coupling to an external physics code or in respose to discrete events.
        for (auto i = 0U; i < n_steps; ++i) {
            // Calculate current time
            const auto t = static_cast<double>(i) * time_step;

            // Apply load in y direction
            interface.turbine.floating_platform.node.loads[1] = 1e6 * sin(2. * M_PI / 20. * t);

            // Apply time varying buoyancy force
            interface.turbine.floating_platform.node.loads[2] =
                buoyancy_force + 0.5 * initial_spring_force * sin(2. * M_PI / 20. * t);

            // Apply time varying moments to platform node
            interface.turbine.floating_platform.node.loads[3] =
                5.0e5 * sin(2. * M_PI / 15. * t);  // rx
            interface.turbine.floating_platform.node.loads[4] =
                1.0e6 * sin(2. * M_PI / 30. * t);  // ry
            interface.turbine.floating_platform.node.loads[5] =
                2.0e7 * sin(2. * M_PI / 60. * t);  // rz

            // Take a single time step
            [[maybe_unused]] const auto converged = interface.Step();
            assert(converged);
        }
    }
    // Make sure to call finalize after all Kynema data structures are deleted
    // and you're ready to exit your application.
    Kokkos::finalize();
    return 0;
}
