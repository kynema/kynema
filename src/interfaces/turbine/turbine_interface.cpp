#include "turbine_interface.hpp"

#include <filesystem>
#include <numbers>
#include <string>

#include "interfaces/components/solution_input.hpp"
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
                controller_input.input_file_path, controller_input.simulation_name,
                controller_input.yaw_control_enabled
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
            outputs_config.output_state_prefixes, outputs_config.buffer_size
        );

        // Write initial state
        this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

        // Write initial time-series data (test values)
        this->WriteTimeSeriesData();
    }
}

void TurbineInterface::UpdateAerodynamicLoads(
    double fluid_density,
    const std::function<std::array<double, 3>(const std::array<double, 3>&)>& inflow_function
) {
    // Get the inflow velocity at the hub node
    this->hub_inflow = inflow_function(
        {this->turbine.hub_node.position[0], this->turbine.hub_node.position[1],
         this->turbine.hub_node.position[2]}
    );

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

    using namespace std::string_literals;
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

    // Yaw angle
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "YawPzn (deg)", this->turbine.yaw_control * 180. / std::numbers::pi
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
        const auto position = this->turbine.azimuth_node.position;
        const auto q_global_local =
            Eigen::Quaternion<double>(position[3], position[4], position[5], position[6]).inverse();
        const auto shaft_loads =
            Eigen::Matrix<double, 3, 1>(this->turbine.shaft_base_to_azimuth.loads.data());
        const auto shaft_forces = q_global_local._transformVector(shaft_loads);
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "RotThrust (kN)", shaft_forces[0] / 1000.
        );
    }

    // Blade based data
    for (auto i : std::views::iota(0U, this->turbine.blades.size())) {
        // Get rotation from global to blade local coordinates
        const auto position = this->turbine.blades[i].nodes[0].position;
        const auto rotation = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(-90. * deg_to_rad, Eigen::Matrix<double, 3, 1>::Unit(1))
        );
        const auto orientation =
            Eigen::Quaternion<double>(position[3], position[4], position[5], position[6]).inverse();
        const auto q_global_to_local = rotation * orientation;

        // Blade root forces in blade coordinates
        const auto blade_root_forces = q_global_to_local._transformVector(
            Eigen::Matrix<double, 3, 1>(this->turbine.blade_pitch[i].loads.data())
        );

        // Blade root forces angles
        const auto blade_number = std::to_string(i + 1);
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "RootFxr (N)"s, blade_root_forces[0]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "RootFyr (N)"s, blade_root_forces[1]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "RootFzr (N)"s, blade_root_forces[2]
        );

        // Blade root moments in blade coordinates
        const auto blade_root_moments = q_global_to_local._transformVector(
            Eigen::Matrix<double, 3, 1>(&this->turbine.blade_pitch[i].loads[3])
        );

        // Blade root moments angles
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "RootMxr (N-m)"s, blade_root_moments[0]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "RootMyr (N-m)"s, blade_root_moments[1]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "RootMzr (N-m)"s, blade_root_moments[2]
        );

        // Blade pitch angles
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "BldPitch"s + blade_number + " (deg)"s,
            this->turbine.blade_pitch_control[i] * 180. / std::numbers::pi
        );

        // Blade tip translational velocity in inertial frame
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "TipTVXg (m_s)",
            this->turbine.blades[i].nodes.back().velocity[0]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "TipTVYg (m_s)",
            this->turbine.blades[i].nodes.back().velocity[1]
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "TipTVZg (m_s)",
            this->turbine.blades[i].nodes.back().velocity[2]
        );

        // Blade tip rotational velocity in inertial frame
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "TipRVXg (deg_s)",
            this->turbine.blades[i].nodes.back().velocity[3] / deg_to_rad
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "TipRVYg (deg_s)",
            this->turbine.blades[i].nodes.back().velocity[4] / deg_to_rad
        );
        this->outputs->WriteValueAtTimestep(
            this->state.time_step, "B"s + blade_number + "TipRVZg (deg_s)",
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

    // Hub inflow velocities in inertial frame
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "WindHubVelX (m_s)", this->hub_inflow[0]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "WindHubVelY (m_s)", this->hub_inflow[1]
    );
    this->outputs->WriteValueAtTimestep(
        this->state.time_step, "WindHubVelZ (m_s)", this->hub_inflow[2]
    );

    // Aerodynamic data
    if (this->aerodynamics) {
        const auto n_blades = 0U;  // this->aerodynamics->bodies.size()
        for (auto i : std::views::iota(0U, n_blades)) {
            const auto& body = this->aerodynamics->bodies[i];
            for (auto j : std::views::iota(0U, body.loads.size())) {
                // Construct the node label
                const auto blade_number = std::to_string(i + 1);
                const auto node_number = std::to_string(j + 1);
                const auto extra_zeros = std::string(3 - node_number.size(), '0');
                auto node_label = "AB"s;
                node_label += blade_number;
                node_label += "N"s;
                node_label += extra_zeros;
                node_label += node_number;

                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Vrel (m_s)",
                    Eigen::Matrix<double, 3, 1>(body.v_rel[j].data()).norm()
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Alpha (deg)", body.alpha[j] / deg_to_rad
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Cn (-)", body.cn[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Ct (-)", body.ct[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Cm (-)", body.cm[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Fxi (N_m)",
                    body.loads[j][0] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Fyi (N_m)",
                    body.loads[j][1] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Fzi (N_m)",
                    body.loads[j][2] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Mxi (N_m)",
                    body.loads[j][3] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Myi (N_m)",
                    body.loads[j][4] / body.delta_s[j]
                );
                this->outputs->WriteValueAtTimestep(
                    this->state.time_step, node_label + "Mzi (N_m)",
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
    controller->io.time = 0.;                                    // Current time (seconds)
    controller->io.azimuth_angle = turbine_input.azimuth_angle;  // Initial azimuth
    controller->io.pitch_blade1_actual = this->turbine.blade_pitch_control[0];  // Blade pitch (rad)
    controller->io.pitch_blade2_actual = this->turbine.blade_pitch_control[1];  // Blade pitch (rad)
    controller->io.pitch_blade3_actual = this->turbine.blade_pitch_control[2];  // Blade pitch (rad)
    controller->io.generator_speed_actual =
        turbine_input.rotor_speed * turbine_input.gear_box_ratio;  // Generator speed (rad/s)
    controller->io.generator_torque_actual =
        turbine_input.generator_power /
        (turbine_input.rotor_speed * turbine_input.gear_box_ratio);         // Generator torque
    controller->io.generator_power_actual = turbine_input.generator_power;  // Generator power (W)
    controller->io.rotor_speed_actual = turbine_input.rotor_speed;          // Rotor speed (rad/s)
    controller->io.horizontal_wind_speed = turbine_input.hub_wind_speed;    // Hub wind speed (m/s)
    controller->io.yaw_angle_actual = this->turbine.yaw_control;            // Yaw angle (rad)

    // Signal first call to controller
    controller->io.status = 0;

    // Make first call to controller to initialize
    controller->CallController();

    this->turbine.torque_control = controller->io.generator_torque_command;
    this->turbine.blade_pitch_control[0] = turbine_input.blade_pitch_angle;
    this->turbine.blade_pitch_control[1] = turbine_input.blade_pitch_angle;
    this->turbine.blade_pitch_control[2] = turbine_input.blade_pitch_angle;
}

void TurbineInterface::ApplyController(double t) {
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
    controller->io.horizontal_wind_speed = sqrt(
        (this->hub_inflow[0] * this->hub_inflow[0]) + (this->hub_inflow[1] * this->hub_inflow[1]) +
        (this->hub_inflow[2] * this->hub_inflow[2])
    );
    controller->io.generator_torque_actual = generator_torque;
    controller->io.generator_power_actual = generator_speed * generator_torque;
    controller->io.pitch_blade1_actual = this->turbine.blade_pitch_control[0];
    controller->io.pitch_blade2_actual = this->turbine.blade_pitch_control[1];
    controller->io.pitch_blade3_actual = this->turbine.blade_pitch_control[2];

    // Loop through blades and calculate out of plane root bending moments
    for (auto i : std::views::iota(0U, this->turbine.blades.size())) {
        // Get rotation from global to blade root coordinates
        // Apex node orientation is the same as blade root node without pitch angle
        const auto position = this->turbine.apex_nodes[i].position;
        const auto q_global_to_local =
            Eigen::Quaternion<double>(position[3], position[4], position[5], position[6]).inverse();

        // Rotate blade root moments into blade root coordinates
        const auto blade_root_moments = q_global_to_local._transformVector(
            Eigen::Matrix<double, 3, 1>(&this->turbine.blade_pitch[i].loads[3])
        );

        // Set out-of-plane root bending moment for each blade (y-axis in blade coords)
        if (i == 0) {
            controller->io.out_of_plane_root_bending_moment_blade1 = blade_root_moments[1];
        } else if (i == 1) {
            controller->io.out_of_plane_root_bending_moment_blade2 = blade_root_moments[1];
        } else if (i == 2) {
            controller->io.out_of_plane_root_bending_moment_blade3 = blade_root_moments[1];
        }
    }

    // Call the controller
    controller->CallController();

    this->turbine.torque_control = controller->io.generator_torque_command;
    this->turbine.blade_pitch_control[0] = controller->io.pitch_collective_command;
    this->turbine.blade_pitch_control[1] = controller->io.pitch_collective_command;
    this->turbine.blade_pitch_control[2] = controller->io.pitch_collective_command;
    this->turbine.yaw_control = controller->YawAngleCommand();
}

void TurbineInterface::OpenOutputFile() {
    if (this->outputs) {
        this->outputs->Open();
    }
}

void TurbineInterface::CloseOutputFile() {
    if (this->outputs) {
        this->outputs->Close();
    }
}

void TurbineInterface::WriteOutput() {
    assert(this->outputs);
    // Write outputs and increment timestep counter

    auto output_region = Kokkos::Profiling::ScopedRegion("Output Data");
    // Write node state outputs
    this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

    // Calculate rotor azimuth and speed -> write rotor time-series data
    this->WriteTimeSeriesData();
}
}  // namespace kynema::interfaces
