#include "turbine_interface.hpp"

#include <filesystem>
#include <format>
#include <numbers>

#include "interfaces/components/solution_input.hpp"
#include "math/quaternion_operations.hpp"
#include "state/clone_state.hpp"
#include "state/copy_state_data.hpp"
#include "step/step.hpp"

namespace kynema::interfaces {

TurbineInterface::TurbineInterface(
    const components::SolutionInput& solution_input, const components::TurbineInput& turbine_input,
    const components::AerodynamicsInput& aerodynamics_input,
    const components::ControllerInput& controller_input,
    const components::OutputsConfig& outputs_config
)
    : model(Model(solution_input.gravity)),
      turbine(turbine_input, model),
      state(model.CreateState<DeviceType>()),
      elements(model.CreateElements<DeviceType>()),
      constraints(model.CreateConstraints<DeviceType>()),
      parameters(
          solution_input.dynamic_solve, solution_input.max_iter, solution_input.time_step,
          solution_input.rho_inf, solution_input.absolute_error_tolerance,
          solution_input.relative_error_tolerance
      ),
      solver(CreateSolver(state, elements, constraints)),
      state_save(CloneState(state)),
      host_state(state),
      host_constraints(constraints) {
    if (aerodynamics_input.is_enabled) {
        auto aero_inputs = std::vector<components::AerodynamicBodyInput>{};
        const auto num_turbine_blades = turbine.blades.size();
        for (auto blade : std::views::iota(0U, num_turbine_blades)) {
            auto blade_node_ids = std::vector<size_t>{};
            std::ranges::transform(
                turbine.blades[blade].nodes, std::back_inserter(blade_node_ids),
                [](const auto& node_data) {
                    return node_data.id;
                }
            );
            aero_inputs.emplace_back(
                blade, blade_node_ids,
                aerodynamics_input.aero_inputs[aerodynamics_input.airfoil_map[blade]]
            );
        }
        aerodynamics = std::make_unique<components::Aerodynamics>(aero_inputs, model.GetNodes());
    }
    // Initialize controller if library path is provided
    if (controller_input.IsEnabled()) {
        try {
            controller = std::make_unique<util::TurbineController>(
                controller_input.shared_lib_path, controller_input.function_name,
                controller_input.input_file_path, controller_input.simulation_name
            );

            // Initialize controller with turbine and solution parameters
            InitializeController(turbine_input, solution_input);
        } catch (const std::runtime_error& e) {
            std::cerr << "Warning: Failed to load controller library '"
                      << controller_input.shared_lib_path << "': " << e.what() << "\n";
            std::cerr << "Continuing without controller." << "\n";
        }
    }

    // Update the host state with current node motion
    this->host_state.CopyFromState(this->state);

    // Update the turbine node motion based on the host state
    this->turbine.GetMotion(this->host_state);

    // Initialize NetCDF writer and write mesh connectivity if output is enabled
    if (outputs_config.Enabled()) {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(outputs_config.output_file_path);

        // Write mesh connectivity to YAML file
        model.ExportMeshConnectivityToYAML(
            outputs_config.output_file_path + "/mesh_connectivity.yaml"
        );

        // Initialize outputs with both node state and time-series files
        this->outputs = std::make_unique<Outputs>(
            outputs_config.output_file_path + "/turbine_interface.nc", this->state.num_system_nodes,
            outputs_config.output_file_path + "/turbine_time_series.nc",
            outputs_config.output_state_prefixes, outputs_config.enable_deformation,
            outputs_config.buffer_size
        );

        // Write initial state
        this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

        // Write initial time-series data (test values)
        this->WriteTimeSeriesData();
    }
}

components::Turbine& TurbineInterface::Turbine() {
    return this->turbine;
}

void TurbineInterface::UpdateAerodynamicLoads(
    double fluid_density,
    const std::function<std::array<double, 3>(const std::array<double, 3>&)>& inflow_function
) {
    if (aerodynamics) {
        auto update_region = Kokkos::Profiling::ScopedRegion("Update Aerodynamic Loads");
        aerodynamics->CalculateMotion(host_state);

        aerodynamics->SetInflowFromFunction(inflow_function);

        aerodynamics->CalculateAerodynamicLoads(fluid_density);

        aerodynamics->CalculateNodalLoads();
    }
}

bool TurbineInterface::Step() {
    auto step_resion = Kokkos::Profiling::ScopedRegion("TurbineInterface::Step");
    // Update the host state with current node loads
    {
        auto forces_region = Kokkos::Profiling::ScopedRegion("Update Forces");
        Kokkos::deep_copy(this->host_state.f, 0.);
        this->turbine.SetLoads(this->host_state);
        if (this->aerodynamics) {
            this->aerodynamics->AddNodalLoadsToState(this->host_state);
        }
        this->host_state.CopyForcesToState(this->state);
    }

    // Solve for state at end of step
    auto converged =
        kynema::Step(this->parameters, this->solver, this->elements, this->state, this->constraints);

    // If not converged, return false
    if (!converged) {
        return false;
    }

    {
        auto update_region = Kokkos::Profiling::ScopedRegion("Update Host Values");
        // Update the host state with current node motion
        this->host_state.CopyFromState(this->state);

        // Update the turbine node motion based on the host state
        this->turbine.GetMotion(this->host_state);

        // Update the host constraints with current constraint loads
        this->host_constraints.CopyFromConstraints(this->constraints);

        // Update the turbine constraint loads based on the host constraints
        this->turbine.GetLoads(this->host_constraints);
    }

    // Write outputs and increment timestep counter
    if (this->outputs) {
        auto output_region = Kokkos::Profiling::ScopedRegion("Output Data");
        // Write node state outputs
        this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

        // Calculate rotor azimuth and speed -> write rotor time-series data
        this->WriteTimeSeriesData();
    }

    return true;
}

void TurbineInterface::SaveState() {
    CopyStateData(this->state_save, this->state);
}

void TurbineInterface::RestoreState() {
    // Copy saved state back to current state
    CopyStateData(this->state, this->state_save);

    // Update the host state with current node motion
    this->host_state.CopyFromState(this->state);

    // Update the turbine node motion based on the host state
    this->turbine.GetMotion(this->host_state);
}

void TurbineInterface::WriteTimeSeriesData() const {
    if (!this->outputs) {
        return;
    }

    constexpr auto rpm_to_radps{0.104719755};            // RPM to rad/s
    constexpr auto deg_to_rad{std::numbers::pi / 180.};  // Degrees to radians

    const auto time_step = static_cast<double>(this->state.time_step);
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "Time (s)", time_step * this->parameters.h
    );
    const auto num_iterations = static_cast<double>(this->solver.convergence_err.size());
    this->outputs->WriteValueAtTimestep(this->state.time_step, "ConvIter (-)", num_iterations);
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "ConvError (-)",
        this->solver.convergence_err.empty() ? 0. : this->solver.convergence_err.back()
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "Azimuth (deg)",
        this->CalculateAzimuthAngle() * 180. / std::numbers::pi
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "RotSpeed (rpm)", this->CalculateRotorSpeed() / rpm_to_radps
    );

    // Tower top displacements
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTDxt (m)", this->turbine.yaw_bearing_node.displacement[0]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTDyt (m)", this->turbine.yaw_bearing_node.displacement[1]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTDzt (m)", this->turbine.yaw_bearing_node.displacement[2]
    );

    // Tower top velocities
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTVxp (m_s)", this->turbine.yaw_bearing_node.velocity[0]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTVyp (m_s)", this->turbine.yaw_bearing_node.velocity[1]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTVzp (m_s)", this->turbine.yaw_bearing_node.velocity[2]
    );

    // Tower top accelerations
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTAxp (m_s^2)", this->turbine.yaw_bearing_node.acceleration[0]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTAyp (m_s^2)", this->turbine.yaw_bearing_node.acceleration[1]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawBrTAzp (m_s^2)", this->turbine.yaw_bearing_node.acceleration[2]
    );

    // Tower base forces
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "TwrBsFxt (kN)", this->turbine.tower_base.loads[0] / 1000.
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "TwrBsFyt (kN)", this->turbine.tower_base.loads[1] / 1000.
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "TwrBsFzt (kN)", this->turbine.tower_base.loads[2] / 1000.
    );

    // Tower base moments
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "TwrBsMxt (kN-m)", this->turbine.tower_base.loads[3] / 1000.
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "TwrBsMyt (kN-m)", this->turbine.tower_base.loads[4] / 1000.
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "TwrBsMzt (kN-m)", this->turbine.tower_base.loads[5] / 1000.
    );

    // Rotor thrust
    {
        const auto q_global_local =
            math::QuaternionInverse(std::span{this->turbine.azimuth_node.position}.subspan<3, 4>());
        const auto shaft_forces = math::RotateVectorByQuaternion(
            q_global_local, std::span{this->turbine.shaft_base_to_azimuth.loads}.subspan<0, 3>()
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "RotThrust (kN)", shaft_forces[0] / 1000.
        );
    }

    // Blade based data
    for (auto i : std::views::iota(0U, this->turbine.blades.size())) {
        // Get rotation from global to blade local coordinates
        const auto q_global_to_local = math::QuaternionCompose(
            math::RotationVectorToQuaternion({0., -90. * deg_to_rad, 0.}),
            math::QuaternionInverse(
                std::span{this->turbine.blades[i].nodes[0].position}.subspan<3, 4>()
            )
        );

        // Blade root forces in blade coordinates
        const auto blade_root_forces = math::RotateVectorByQuaternion(
            q_global_to_local, std::span{this->turbine.blade_pitch[i].loads}.subspan<0, 3>()
        );

        // Blade root forces angles
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}RootFxr (N)", i + 1), blade_root_forces[0]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}RootFyr (N)", i + 1), blade_root_forces[1]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}RootFzr (N)", i + 1), blade_root_forces[2]
        );

        // Blade root moments in blade coordinates
        const auto blade_root_moments = math::RotateVectorByQuaternion(
            q_global_to_local, std::span{this->turbine.blade_pitch[i].loads}.subspan<3, 3>()
        );

        // Blade root moments angles
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}RootMxr (N-m)", i + 1), blade_root_moments[0]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}RootMyr (N-m)", i + 1), blade_root_moments[1]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}RootMzr (N-m)", i + 1), blade_root_moments[2]
        );

        // Blade pitch angles
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("BldPitch{} (deg)", i + 1),
            this->turbine.blade_pitch_control[i] * 180. / std::numbers::pi
        );

        // Blade tip translational velocity in inertial frame
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}TipTVXg (m_s)", i + 1),
            this->turbine.blades[i].nodes.back().velocity[0]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}TipTVYg (m_s)", i + 1),
            this->turbine.blades[i].nodes.back().velocity[1]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}TipTVZg (m_s)", i + 1),
            this->turbine.blades[i].nodes.back().velocity[2]
        );

        // Blade tip rotational velocity in inertial frame
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}TipRVXg (deg_s)", i + 1),
            this->turbine.blades[i].nodes.back().velocity[3] / deg_to_rad
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}TipRVYg (deg_s)", i + 1),
            this->turbine.blades[i].nodes.back().velocity[4] / deg_to_rad
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, std::format("B{}TipRVZg (deg_s)", i + 1),
            this->turbine.blades[i].nodes.back().velocity[5] / deg_to_rad
        );
    }

    // Generator torque and power if controller is present
    if (this->controller) {
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "GenTq (kN-m)",
            this->controller->io.generator_torque_command / 1000.
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "GenPwr (kW)", this->controller->io.generator_power_actual / 1000.
        );
    }

    // Aerodynamic data
    if (this->aerodynamics) {
        for (auto i : std::views::iota(0U, this->aerodynamics->bodies.size())) {
            const auto& body = this->aerodynamics->bodies[i];
            for (auto j : std::views::iota(0U, body.loads.size())) {
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Vrel (m_s)", i + 1, j + 1),
                    math::Norm(body.v_rel[j])
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Alpha (deg)", i + 1, j + 1),
                    body.alpha[j] / deg_to_rad
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Cn (-)", i + 1, j + 1), body.cn[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Ct (-)", i + 1, j + 1), body.ct[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Cm (-)", i + 1, j + 1), body.cm[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Fxi (N_m)", i + 1, j + 1),
                    body.loads[j][0] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Fyi (N_m)", i + 1, j + 1),
                    body.loads[j][1] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Fzi (N_m)", i + 1, j + 1),
                    body.loads[j][2] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Mxi (N_m)", i + 1, j + 1),
                    body.loads[j][3] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Myi (N_m)", i + 1, j + 1),
                    body.loads[j][4] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, std::format("AB{}N{:03}Mzi (N_m)", i + 1, j + 1),
                    body.loads[j][5] / body.delta_s[j]
                );
            }
        }
    }
}

double TurbineInterface::CalculateAzimuthAngle() const {
    const auto azimuth_constraint_id = this->turbine.shaft_base_to_azimuth.id;
    double azimuth = this->constraints.host_output(azimuth_constraint_id, 0);

    // Normalize azimuth angle to range [0, 2π) radians
    azimuth = std::fmod(azimuth, 2. * std::numbers::pi);
    if (azimuth < 0) {
        azimuth += 2. * std::numbers::pi;
    }

    return azimuth;
}

double TurbineInterface::CalculateRotorSpeed() const {
    const auto azimuth_constraint_id = this->turbine.shaft_base_to_azimuth.id;
    return this->constraints.host_output(azimuth_constraint_id, 1);
}

std::array<double, 3> TurbineInterface::GetHubNodePosition() const {
    return std::array{
        turbine.hub_node.position[0], turbine.hub_node.position[1], turbine.hub_node.position[2]
    };
}

void TurbineInterface::InitializeController(
    const components::TurbineInput& turbine_input, const components::SolutionInput& solution_input
) {
    if (!controller) {
        return;
    }

    // Set controller constant parameters
    controller->io.dt = solution_input.time_step;           // Time step size (seconds)
    controller->io.pitch_actuator_type_req = 0;             // Pitch position actuator
    controller->io.pitch_control_type = 0;                  // Collective pitch control
    controller->io.n_blades = turbine_input.blades.size();  // Number of blades

    // Set controller initial values
    controller->io.time = 0.;                                              // Current time (seconds)
    controller->io.azimuth_angle = turbine_input.azimuth_angle;            // Initial azimuth
    controller->io.pitch_blade1_actual = turbine_input.blade_pitch_angle;  // Blade pitch (rad)
    controller->io.pitch_blade2_actual = turbine_input.blade_pitch_angle;  // Blade pitch (rad)
    controller->io.pitch_blade3_actual = turbine_input.blade_pitch_angle;  // Blade pitch (rad)
    controller->io.generator_speed_actual =
        turbine_input.rotor_speed * turbine_input.gear_box_ratio;  // Generator speed (rad/s)
    controller->io.generator_torque_actual =
        turbine_input.generator_power /
        (turbine_input.rotor_speed * turbine_input.gear_box_ratio);         // Generator torque
    controller->io.generator_power_actual = turbine_input.generator_power;  // Generator power (W)
    controller->io.rotor_speed_actual = turbine_input.rotor_speed;          // Rotor speed (rad/s)
    controller->io.horizontal_wind_speed = turbine_input.hub_wind_speed;    // Hub wind speed (m/s)

    // Signal first call to controller
    controller->io.status = 0;

    // Make first call to controller to initialize
    controller->CallController();

    this->turbine.torque_control = controller->io.generator_torque_command;
    this->turbine.blade_pitch_control[0] = turbine_input.blade_pitch_angle;
    this->turbine.blade_pitch_control[1] = turbine_input.blade_pitch_angle;
    this->turbine.blade_pitch_control[2] = turbine_input.blade_pitch_angle;
}

void TurbineInterface::ApplyController(double t, double hub_wind_speed) {
    if (!controller) {
        return;
    }

    // Update controller inputs from current system state
    // Update time and azimuth
    controller->io.status = 1;
    controller->io.time = t;
    controller->io.azimuth_angle = CalculateAzimuthAngle();

    // Update rotor and generator speeds
    const double rotor_speed = CalculateRotorSpeed();
    controller->io.rotor_speed_actual = rotor_speed;
    controller->io.generator_speed_actual =
        rotor_speed * this->turbine.GetTurbineInput().gear_box_ratio;

    // Update generator power and torque
    const double generator_speed = controller->io.generator_speed_actual;
    const double generator_torque = this->turbine.torque_control;
    controller->io.horizontal_wind_speed = hub_wind_speed;
    controller->io.generator_torque_actual = generator_torque;
    controller->io.generator_power_actual = generator_speed * generator_torque;
    controller->io.pitch_blade1_actual = this->turbine.blade_pitch_control[0];
    controller->io.pitch_blade2_actual = this->turbine.blade_pitch_control[1];
    controller->io.pitch_blade3_actual = this->turbine.blade_pitch_control[2];

    // Call the controller
    controller->CallController();

    this->turbine.torque_control = controller->io.generator_torque_command;
    this->turbine.blade_pitch_control[0] = controller->io.pitch_collective_command;
    this->turbine.blade_pitch_control[1] = controller->io.pitch_collective_command;
    this->turbine.blade_pitch_control[2] = controller->io.pitch_collective_command;
}
}  // namespace kynema::interfaces
