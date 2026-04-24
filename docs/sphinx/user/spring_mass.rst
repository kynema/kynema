Example: Spring-Mass System
===========================

This example will walk through how to run a simulation of a chain of masses linked together by linear springs and anchored at either end.
We'll use Kynema-FMB's low level API, which, unlike Kynema-FMB's high level APIs, will require us to manually set up all nodes and their connectivities.
This extra complexity is the trade-off required for unlimited freedom.
For the most up to date and working version of this code, see ``tests/documentation_tests/spring_mass_system/``.

As with any C++ program, start with the includes.
As a Kokkos-based library, you'll need to include ``Kokkos_Core.hpp`` for setup, teardown, and working with Kynema-FMB's data structures.
From Kynema-FMB, you'll have to include ``model.hpp`` for the Model class, our tool for setting up and creating the system, and ``step.hpp`` for the Step function which performs the action of system assembly and solve.

.. code-block:: cpp

    #include <array>
    #include <cassert>
    #include <Kokkos_Core.hpp>
    #include <model/model.hpp>
    #include <step/step.hpp>

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

Now, we define the mass matrix.

.. code-block:: cpp

    constexpr auto mass = 1.;
    constexpr auto inertia = 1.;
    constexpr auto mass_matrix = std::array{
        std::array{mass, 0.,   0.,   0.,      0.,      0.     },
        std::array{0.,   mass, 0.,   0.,      0.,      0.     },
        std::array{0.,   0.,   mass, 0.,      0.,      0.     },
        std::array{0.,   0.,   0.,   inertia, 0.,      0.     },
        std::array{0.,   0.,   0.,   0.,      inertia, 0.     },
        std::array{0.,   0.,   0.,   0.,      0.,      inertia},
    };

A Model is Kynema-FMB's low level interface for specifying elements, nodes, constraints, and their connectivities.
One everything has been specified, we will use model to create Kynema-FMB's fundamental data structures and advance the problem in time.

.. code-block:: cpp

    auto model = kynema_fmb::Model();

To add a node, we call the AddNode method on Model, which creates a NodeBuilder object.
This factory lets us string together function calls to specify the initial position, velocity, and acceleration in a human readable fashion.
Once we finalized by calling the Build method, the NodeBuilder adds a node to the model and returns its newly created ID number.
This ID will be used for specifying elements and constraints.

For this problem, we'll add a series of equally spaced nodes for each mass and an anchor point at the beginning and end of the list.
We'll store the node-ids for each of these for later use defining the physics and connectivity information.

.. code-block:: cpp

    constexpr auto number_of_masses = 10U;
    constexpr auto displacement = 0.5;
    auto anchor_node_ids = std::array<size_t, 2>{};
    auto mass_node_ids = std::array<size_t, number_of_masses>{};
    auto position = 0.;
    anchor_node_ids.front() = model.AddNode().SetPosition(position, 0., 0., 1., 0., 0., 0.).Build();
    for (auto& mass_node_id : mass_node_ids) {
        position += displacement;
        mass_node_id = model.AddNode().SetPosition(position, 0., 0., 1., 0., 0., 0.).Build();
    }
    position += displacement;
    anchor_node_ids.back() = model.AddNode().SetPosition(position, 0., 0., 1., 0., 0., 0.).Build();

To add a mass element to the model, we will need a single node id number and the mass matrix containing the mass and inertia information.

We'll add a mass element for each of the mass nodes

.. code-block:: cpp

    for (auto mass_node_id : mass_node_ids) {
        model.AddMassElement(mass_node_id, mass_matrix);
    }

To add a spring to the model, we will need a node id number, a stiffness, and an undisplaced length where the spring force is zero.

We'll add a spring element between each of the mass elements and between the ending mass elements and our anchor points.

.. code-block:: cpp

    const auto stiffness = 10.;
    const auto length = 0.;
    model.AddSpringElement(anchor_node_ids.front(), mass_node_ids.front(), stiffness, length);
    for (auto index = 0U; index < number_of_masses-1; ++index) {
        model.AddSpringElement(mass_node_ids[index], mass_node_ids[index+1U], stiffness, length);
    }
    model.AddSpringElement(mass_node_ids.back(), anchor_node_ids.back(), stiffness, length);

Each of the anchor nodes requires a fixed boundary condition, which will prevent it from either moving or rotating.

.. code-block:: cpp

    for (auto anchor_node_id : anchor_node_ids) {
        model.AddFixedBC(anchor_node_id);
    }

Now that the problem has been fully described in the model, we will create Kynema-FMB's main data structures: State, Elements, Constraints, and Solver.

The CreateSystem method takes an optional template argument with a Kokkos device describing where the system will reside and run.
By default, it uses Kokkos' default execution/memory space, so a serial build will run on the CPU, a CUDA build will run on a CUDA device, etc.

The CreateSolver<> function uses the connectivity defined in the State, Elements, and Constraints structures to construct the Solver object.

State contains the current state (position, velocity, etc) information for each node.

Elements contains each a Beams, Masses, and Springs structure.
These contain the connectivity and basis information or all of the elements of the respective type.

Constraints contains the connectivity information for each constraint in the system.

Solver contains the linear system (sparse matrix, RHS) and linear system solver

.. code-block:: cpp

    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = kynema_fmb::CreateSolver<>(state, elements, constraints);

The final stage is to create a StepParameters object, which contains information like the number of non-linear iterations, time step size, and numerical damping factor used to take a single time step.

.. code-block:: cpp

    constexpr auto num_steps = 1000;
    constexpr bool is_dynamic_solve(true);
    constexpr size_t max_iter(6);
    constexpr double rho_inf(0.);
    const double final_time = 2. * M_PI * sqrt(mass / stiffness);
    const double step_size(final_time / static_cast<double>(num_steps));
    auto parameters = kynema_fmb::StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

Kynema-FMB allows the user to control the actual time stepping process.
This includes setting forces, post-processing data, or coupling to other codes.
In this example, we'll check that none of the nodes have moved - the chain is in constant tension and equilibrium.

The current state is stored in the State object's q member.
This is a Kokkos view of size num_nodes x 7.
This view lives on device, so we can't access it directly from host code.
Here, we create a mirror view on host and, at each time step, copy the data to host and check the value of x-displacement at each node.

.. code-block:: cpp

    auto q = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.q);
    for (auto i = 0; i < 400; ++i) {
        [[maybe_unused]] const auto converged =
            kynema_fmb::Step(parameters, solver, elements, state, constraints);
        assert(converged);
        Kokkos::deep_copy(q, state.q);
        for (auto node = 0; node < 7; ++node) {
            assert(std::abs(q(node, 0)) < 1.e-14);
        }
    }
