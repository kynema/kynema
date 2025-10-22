Example: Heavy Top Problem
==========================

This example will walk through how to run a simulation of a processing top using Kynema's low level API.
Unline Kynema's high level APIs, you will have to manually set up all nodes and their connectivities.
While this extra work adds complexity compared to the higher level APIs, it also provides unlimited freedom.
The heavy top problem is one of the simplest problems you'll want to solve with Kynema, so it is a good introduction to our low level APIs.
For the most up to date and working version of this code, see ``tests/documentation_tests/heavy_top/``.

As with any C++ program, start with the includes.
To set up problems like this one, you might need to perform some linear algebra - we'll include ``Eigen/Dense`` to perform these operations, but you can use any libraries or hand written code as you see fit.
As a Kokkos-based library, you'll need to include ``Kokkos_Core.hpp`` for setup, teardown, and working with Kynema's data structures.
From Kynema, you'll have to include ``model.hpp`` for the Model class, our tool for setting up and creating the system, and ``step.hpp`` for the Step function which performs the action of system asembly and solve.

.. code-block:: cpp

    #include <array>
    #include <cassert>

    #include <Eigen/Dense>
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

Now, we define the mass matrix and initial position, velocity, and acceleration.

.. code-block:: cpp

    constexpr auto mass = 15.;
    constexpr auto inertia = std::array{0.234375, 0.46875, 0.234375};
    const auto x = Eigen::Matrix<double, 3, 1>(0., 1., 0.);
    const auto omega = Eigen::Matrix<double, 3, 1>(0., 150., -4.61538);
    const auto x_dot = omega.cross(x);
    const auto omega_dot = std::array{661.3461692307691919, 0., 0.};
    const auto x_ddot = std::array{0., -21.3017325444000001, -30.9608307692308244};

A Model is Kynema's low level interface for specifying elements, nodes, constraints, and their connectivities.
One everything has been specified, we will use model to create Kynema's fundamental data structures and advance the problem in time.

.. code-block:: cpp

    auto model = kynema::Model();

To add a node, we call the AddNode method on Model, which creates a NodeBuilder object.
This factory lets us string together function calls to specify the initial position, velocity, and acceleration in a human readable fashion.
Once we finalized by calling the Build method, the NodeBuilder adds a node to the model and returns its newly created ID number.
This ID will be used for specifying elements and constraints.

For this problem, we'll add two nodes: one for the mass and one at the origin for defining constraints.

.. code-block:: cpp

    auto mass_node_id =
        model.AddNode()
            .SetPosition(x[0], x[1], x[2], 1., 0., 0., 0.)
            .SetVelocity(x_dot[0], x_dot[1], x_dot[2], omega[0], omega[1], omega[2])
            .SetAcceleration(x_ddot[0], x_ddot[1], x_ddot[2], omega_dot[0], omega_dot[1], omega_dot[2])
            .Build();
    auto ground_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();

To add a mass element to the model, we will need a single node id number and the mass matrix containing the mass and inertia information.

.. code-block:: cpp

    model.AddMassElement(
        mass_node_id, {{
                          {mass, 0., 0., 0., 0., 0.},
                          {0., mass, 0., 0., 0., 0.},
                          {0., 0., mass, 0., 0., 0.},
                          {0., 0., 0., inertia[0], 0., 0.},
                          {0., 0., 0., 0., inertia[1], 0.},
                          {0., 0., 0., 0., 0., inertia[2]},
                      }}
    );

This problem requires two constraints: a rigid joint prescribing that the center of mass remains a constant distance from the ground node and a prescribed boundary condition forcing the ground node to remain stationary.

.. code-block:: cpp

    model.AddRigidJoint6DOFsTo3DOFs(std::array{mass_node_id, ground_node_id});
    model.AddPrescribedBC3DOFs(ground_node_id);

The gravity vector for the problem is set using the well named SetGravity method

.. code-block:: cpp

    model.SetGravity(0., 0., -9.81);

Now that the problem has been fully described in the model, we will create Kynema's main data structures: State, Elements, Constraints, and Solver.
The CreateSystemWithSolver<> method takes an optional template argument with a Kokkos device describing where the system will reside and run.
By default, it uses Kokkos' default execution/memory space, so a serial build will run on the CPU, a CUDA build will run on a CUDA device, etc.

State contains the current state (position, velocity, etc) information for each node.

Elements contains each a Beams, Masses, and Springs structure.
These contain the connectivity and basis information or all of the elements of the respective type.

Constraints contains the connectivity information for each constraint in the system.

Solver contains the linear system (sparse matrix, RHS) and linear system solver

.. code-block:: cpp

    auto [state, elements, constraints, solver] = model.CreateSystemWithSolver<>();

The final stage is to create a StepParameters object, which contains information like the number of non-linear iterations, time step size, and numerical damping factor used to take a single time step.

.. code-block:: cpp

    constexpr auto is_dynamic_solve(true);
    constexpr auto max_iter(10UL);
    constexpr auto step_size(0.002);
    constexpr auto rho_inf(0.9);
    constexpr auto a_tol(1e-5);
    constexpr auto r_tol(1e-3);
    auto parameters = kynema::StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf, a_tol, r_tol);

Kynema allows the user to control the actual time stepping process.
This includes setting forces, post-processing data, coupling to other codes.
This example does none of that.
At each time step, we call Kynema's Step function and pass in the previously created structures.

.. code-block:: cpp

    for (auto i = 0; i < 400; ++i) {
        const auto converged = kynema::Step(parameters, solver, elements, state, constraints);
    }

Finally, we can check that our solution is correct.
The current state is stored in the State object's q member.
This is a Kokkos view of size num_nodes x 7.
This View lives on device, so we can't access it directly from host code.
Here, we create a mirror view on host and then check the values.
For more information on working with Kokkos data structures, see the Kokkos documentation.

.. code-block:: cpp

    const auto q = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, state.q);
    assert(std::abs(q(0, 0) - -0.42217802273894345) < 1e-10);
    assert(std::abs(q(0, 1) - -0.09458263530050703) < 1e-10);
    assert(std::abs(q(0, 2) - -0.04455460488952848) < 1e-10);
    assert(std::abs(q(0, 3) - -0.17919607435565366) < 1e-10);
    assert(std::abs(q(0, 4) - 0.21677896640311572) < 1e-10);
    assert(std::abs(q(0, 5) - -0.95947769608535960) < 1e-10);
    assert(std::abs(q(0, 6) - -0.017268392381761217) < 1e-10);
