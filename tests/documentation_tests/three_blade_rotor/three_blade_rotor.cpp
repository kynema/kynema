#include <array>
#include <cassert>

#include <Eigen/Dense>
#include <Kokkos_Core.hpp>
#include <model/model.hpp>
#include <step/step.hpp>

int main() {
    // Kynema is based on Kokkos for performance portability.  Make sure to
    // call Kokkos::initialize before creating any Kynema data structures
    // and Kokkos::finalize after all of those data structures have been destroyed.
    Kokkos::initialize();
    {
        // We first define our Beam Sections, which consist of the mass and stiffness
        // matrices defined at given nodes along the beam.  These physical properties
        // will be interpolated to the quadrature points for use in our simulation.
        // This problem uses a constant-cross-section blade, so the sections at the root
        // and tip of the beam will be identical, but any number of complex cross-sections
        // can be modeled
        constexpr auto mass_matrix = std::array{
            std::array{8.538e-2, 0., 0., 0., 0., 0.},   std::array{0., 8.538e-2, 0., 0., 0., 0.},
            std::array{0., 0., 8.538e-2, 0., 0., 0.},   std::array{0., 0., 0., 1.4433e-2, 0., 0.},
            std::array{0., 0., 0., 0., 0.40972e-2, 0.}, std::array{0., 0., 0., 0., 0., 1.0336e-2},
        };
        constexpr auto stiffness_matrix = std::array{
            std::array{1368.17e3, 0., 0., 0., 0., 0.},
            std::array{0., 88.56e3, 0., 0., 0., 0.},
            std::array{0., 0., 38.78e3, 0., 0., 0.},
            std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
            std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
            std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
        };
        const auto sections = std::vector{
            kynema::BeamSection(0., mass_matrix, stiffness_matrix),
            kynema::BeamSection(1., mass_matrix, stiffness_matrix),
        };

        // We now define the node locations where our solution will be defined.  In this case, we
        // will use the a six node GLL basis.  Like the sections, these will be defined along the
        // blade's reference axis.
        const auto node_s = std::array{
            0., 0.11747233803526763, 0.35738424175967748, 0.64261575824032247, 0.88252766196473242,
            1.
        };

        // As a final piece of information, we will define the quadrature points and their weights.
        // We've precalculated these to be a seven point GL quadrature rule, which will be
        // sufficiently accurate for both our basis and mass/stiffness matrix distribution.
        const auto quadrature = std::vector<std::array<double, 2>>{
            {-0.9491079123427585, 0.1294849661688697},  {-0.7415311855993943, 0.27970539148927664},
            {-0.40584515137739696, 0.3818300505051189}, {6.123233995736766e-17, 0.4179591836734694},
            {0.4058451513773971, 0.3818300505051189},   {0.7415311855993945, 0.27970539148927664},
            {0.9491079123427585, 0.1294849661688697},
        };

        // A Model is Kynema's low level interface for specifying elements, nodes, constraints,
        // and their connectivities.  Once everything has be specified, we will use to model to
        // create Kynema's fundamental data structures and advance the problem in time.
        auto model = kynema::Model();

        // The aptly named SetGravity method is used to set the gravity vector for the problem.
        model.SetGravity(0., 0., -9.81);

        // When specifying the beam elements, we'll also set the initial velocity.  To help
        // formulate this, we specify the rotor velocity (both translational and rotational)
        // and the origin about which we'll rotate.
        constexpr auto velocity = std::array{0., 0., 0., 0., 0., 1.};
        constexpr auto origin = std::array{0., 0., 0.};
        constexpr auto hub_radius = 2.;

        // We'll now define three beam elements to be our main rotor.  Each of these beams will be
        // identical, but we'll rotate each of them by 120 degrees around the origin to create a
        // rotor like one would see on a wind turbine.
        constexpr auto num_blades = 3;
        for (auto blade_number = 0; blade_number < num_blades; ++blade_number) {
            auto beam_node_ids = std::vector<size_t>(node_s.size());
            std::transform(
                std::cbegin(node_s), std::cend(node_s), std::begin(beam_node_ids),
                [&](auto s) {
                    return model.AddNode()
                        .SetElemLocation(s)
                        .SetPosition(10. * s, 0., 0., 1., 0., 0., 0.)
                        .Build();
                }
            );
            auto blade_elem_id = model.AddBeamElement(beam_node_ids, sections, quadrature);
            const auto rotation_angle = 2. * std::numbers::pi * blade_number / num_blades;
            const auto rotation_quaternion = Eigen::Quaternion<double>(
                Eigen::AngleAxis(rotation_angle, Eigen::Matrix<double, 3, 1>::Unit(2))
            );
            const auto rotation_quaternion_array = std::array{
                rotation_quaternion.w(), rotation_quaternion.x(), rotation_quaternion.y(),
                rotation_quaternion.z()
            };
            model.TranslateBeam(blade_elem_id, std::array{hub_radius, 0., 0.});
            model.RotateBeamAboutPoint(blade_elem_id, rotation_quaternion_array, origin);
            model.SetBeamVelocityAboutPoint(blade_elem_id, velocity, origin);
        }

        // To control the rotation of the turbine, we create a node to act as a hub and attach
        // the nearest node of each beam element to the hub with a rigid joint constraint.
        // We then create a prescribed boundary condidition constraint on the hub, which we will
        // modify during time stepping to create rotation.
        auto hub_node_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
        for (const auto& beam_element : model.GetBeamElements()) {
            model.AddRigidJointConstraint(std::array{hub_node_id, beam_element.node_ids.front()});
        }
        auto hub_bc_id = model.AddPrescribedBC(hub_node_id);

        // Now that the problem has been fully described in the model, we will create Kynema's
        // main data structures: State, Elements, Constraints, and Solver.
        //
        // The CreateSystem method takes an optional template argument with a Kokkos device
        // describing where the system will reside and run.  By default, it uses Kokkos' default
        // execution/memory space, so a serial build will run on the CPU, a CUDA build will run on a
        // CUDA device, etc.
        //
        // The CreateSolver<> function uses the connectivity defined in the State, Elements, and
        // Constraints structures to construct the Solver object.
        //
        // State contains the current state (position, velocity, etc) information for each node.
        //
        // Elements contains each a Beams, Masses, and Springs structure.  These contain the
        // connectivity and basis information or all of the elements of the respective type.
        //
        // Constraints contains the connectivity information for each constraint in the system.
        //
        // Solver contains the linear system (sparse matrix, RHS) and linear system solver
        auto [state, elements, constraints] = model.CreateSystem();
        auto solver = kynema::CreateSolver<>(state, elements, constraints);

        // The final stage is to create a StepParameters object, which contains information like
        // the number of non-linear iterations, time step size, and numerical damping factor used
        // to take a single time step.
        const bool is_dynamic_solve(true);
        const int max_iter(4);
        const double step_size(0.01);
        const double rho_inf(0.9);
        const double t_end(0.1);
        const auto num_steps = static_cast<size_t>(std::floor(t_end / step_size + 1.0));
        auto parameters = kynema::StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

        // Kynema allows the user to control the actual time stepping process.  This includes
        // setting forces, post-processing data, or coupling to other codes.
        // For this problem, we will prescribe a rotation on the hub boundary condition, which will
        // be transmitted to the blades through their respective constraints.
        for (auto i = 0U; i < num_steps; ++i) {
            const auto omega_vector = Eigen::Matrix<double, 3, 1>(&velocity[3]);
            const auto omega_axis = omega_vector.normalized();
            const auto omega = omega_vector.norm() * step_size * static_cast<double>(i + 1);
            const auto q_hub =
                Eigen::Quaternion<double>(Eigen::AngleAxis<double>(omega, omega_axis));
            const auto u_hub = std::array{0., 0., 0., q_hub.w(), q_hub.x(), q_hub.y(), q_hub.z()};
            constraints.UpdateDisplacement(hub_bc_id, u_hub);
            [[maybe_unused]] const auto converged =
                kynema::Step(parameters, solver, elements, state, constraints);
            assert(converged);
        }
    }
    // Make sure to call finalize after all Kynema data structures are deleted
    // and you're ready to exit your application.
    Kokkos::finalize();
    return 0;
}
