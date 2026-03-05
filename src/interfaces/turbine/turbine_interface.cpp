#include "turbine_interface.hpp"

#include <algorithm>
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
        // If there are more entries in the airfoil inputs than blades in the turbine
        // Assume that the last mapping index is for the tower
        if (aerodynamics_input.airfoil_map.size() > num_turbine_blades) {
            auto tower_node_ids = std::vector<size_t>{};
            std::ranges::transform(
                turbine.tower.nodes, std::back_inserter(tower_node_ids),
                [](const auto& node_data) {
                    return node_data.id;
                }
            );
            aero_inputs.emplace_back(
                num_turbine_blades, tower_node_ids,
                aerodynamics_input.aero_inputs[aerodynamics_input.airfoil_map.back()]
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
                turbine_input.nacelle_yaw_angle, controller_input.yaw_control_enabled
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

        // Build time-series schema
        this->BuildTimeSeriesSchema();

        // Initialize outputs with the NEW constructor that supports channels
        this->outputs = std::make_unique<Outputs>(
            outputs_config.output_file_path + "/turbine_interface.nc", this->state.num_system_nodes,
            outputs_config.output_file_path + "/turbine_time_series.nc",
            this->time_series_channels_,  // Channel names from schema
            std::vector<std::string>{},   // Empty units (units are in channel names)
            outputs_config.output_state_prefixes,
            outputs_config.buffer_size,  // Node state buffer
            outputs_config.buffer_size   // Time-series buffer
        );

        // Write initial state
        this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

        // Write initial time-series data (test values)
        this->WriteTimeSeriesData();
    }
}

void TurbineInterface::BuildTimeSeriesSchema() {
    //------------------------------------------
    // lambdas
    //------------------------------------------
    // lambda to add a channel to the time-series schema
    auto add = [&](const std::string& name) -> size_t {
        this->time_series_channels_.push_back(name);
        return this->time_series_channels_.size() - 1;
    };

    // lambda to pad a number with leading zeros to 3 digits (e.g., "AB1N001")
    auto pad_3_digits = [](size_t value) -> std::string {
        auto s = std::to_string(value);
        if (s.size() >= 3) {
            return s;  // Already 3 digits or more
        }
        return std::string(3 - s.size(), '0') + s;  // Pad with leading zeros
    };

    //------------------------------------------
    // build the time-series schema
    //------------------------------------------

    this->time_series_channels_.clear();      // clear the channel names
    this->time_series_enabled_ = false;       // disable time-series output
    this->index_map_ = TimeSeriesIndexMap{};  // clear the index map

    // Basic simulation parameter and other misc. channels
    this->index_map_.time_seconds = add("Time (s)");
    this->index_map_.num_convergence_iterations = add("ConvIter (-)");
    this->index_map_.convergence_error = add("ConvError (-)");
    this->index_map_.azimuth_angle_degrees = add("Azimuth (deg)");
    this->index_map_.rotor_speed_rpm = add("RotSpeed (rpm)");
    this->index_map_.yaw_position_degrees = add("YawPzn (deg)");

    // Tower top channels
    this->index_map_.tower_top_displacement_start = add("YawBrTDxt (m)");
    add("YawBrTDyt (m)");
    add("YawBrTDzt (m)");

    this->index_map_.tower_top_velocity_start = add("YawBrTVxp (m_s)");
    add("YawBrTVyp (m_s)");
    add("YawBrTVzp (m_s)");

    this->index_map_.tower_top_acceleration_start = add("YawBrTAxp (m_s^2)");
    add("YawBrTAyp (m_s^2)");
    add("YawBrTAzp (m_s^2)");

    // Tower base loads
    this->index_map_.tower_base_force_start = add("TwrBsFxt (kN)");
    add("TwrBsFyt (kN)");
    add("TwrBsFzt (kN)");

    this->index_map_.tower_base_moment_start = add("TwrBsMxt (kN-m)");
    add("TwrBsMyt (kN-m)");
    add("TwrBsMzt (kN-m)");

    // Rotor thrust
    this->index_map_.rotor_thrust_kN = add("RotThrust (kN)");

    // Blade channels
    this->index_map_.blade_channel_offsets.clear();
    this->index_map_.blade_channel_offsets.reserve(this->turbine.blades.size());
    for (size_t i = 0; i < this->turbine.blades.size(); ++i) {
        const auto blade_number = std::to_string(i + 1);
        this->index_map_.blade_channel_offsets.push_back(this->time_series_channels_.size());

        add("B" + blade_number + "RootFxr (N)");
        add("B" + blade_number + "RootFyr (N)");
        add("B" + blade_number + "RootFzr (N)");
        add("B" + blade_number + "RootMxr (N-m)");
        add("B" + blade_number + "RootMyr (N-m)");
        add("B" + blade_number + "RootMzr (N-m)");
        add("BldPitch" + blade_number + " (deg)");
        add("B" + blade_number + "TipTVXg (m_s)");
        add("B" + blade_number + "TipTVYg (m_s)");
        add("B" + blade_number + "TipTVZg (m_s)");
        add("B" + blade_number + "TipRVXg (deg_s)");
        add("B" + blade_number + "TipRVYg (deg_s)");
        add("B" + blade_number + "TipRVZg (deg_s)");
    }

    // Controller channels
    this->index_map_.has_controller_channels = (this->controller != nullptr);
    if (this->index_map_.has_controller_channels) {
        this->index_map_.generator_torque_kNm = add("GenTq (kN-m)");
        this->index_map_.generator_power_kW = add("GenPwr (kW)");
    }

    // Hub inflow velocities channels
    this->index_map_.hub_inflow_start = add("WindHubVelX (m_s)");
    add("WindHubVelY (m_s)");
    add("WindHubVelZ (m_s)");

    // Aerodynamic channels
    this->index_map_.aero_body_offsets.clear();
    this->index_map_.aero_section_counts.clear();

    if (this->aerodynamics) {
        const size_t n_bodies =
            std::min(this->aerodynamics->bodies.size(), this->turbine.blades.size());

        this->index_map_.aero_body_offsets.reserve(n_bodies);
        this->index_map_.aero_section_counts.reserve(n_bodies);

        for (size_t body_index = 0; body_index < n_bodies; ++body_index) {
            const auto& body = this->aerodynamics->bodies[body_index];

            this->index_map_.aero_body_offsets.push_back(this->time_series_channels_.size());
            this->index_map_.aero_section_counts.push_back(body.loads.size());

            for (size_t section_index = 0; section_index < body.loads.size(); ++section_index) {
                const auto blade_number = std::to_string(body_index + 1);
                const auto node_label = "AB" + blade_number + "N" + pad_3_digits(section_index + 1);

                add(node_label + "Vrel (m_s)");
                add(node_label + "Alpha (deg)");
                add(node_label + "Cn (-)");
                add(node_label + "Ct (-)");
                add(node_label + "Cm (-)");
                add(node_label + "Fxi (N_m)");
                add(node_label + "Fyi (N_m)");
                add(node_label + "Fzi (N_m)");
                add(node_label + "Mxi (N_m)");
                add(node_label + "Myi (N_m)");
                add(node_label + "Mzi (N_m)");
            }
        }
    }

    // Allocate row buffer
    this->time_series_row_buffer_.assign(this->time_series_channels_.size(), 0.);
    this->time_series_enabled_ = true;
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

void TurbineInterface::WriteTimeSeriesData() {
    if (!this->outputs || !this->time_series_enabled_) {
        return;
    }

    // initialize the time-series row buffer
    std::fill(this->time_series_row_buffer_.begin(), this->time_series_row_buffer_.end(), 0.);

    //------------------------------------------
    // conversion constants
    //------------------------------------------
    constexpr auto rpm_to_radps{0.104719755};
    constexpr auto rad_to_deg{180. / std::numbers::pi};
    constexpr auto deg_to_rad{std::numbers::pi / 180.};

    //------------------------------------------
    // write the time-series data
    //------------------------------------------

    const size_t time_step{this->state.time_step};

    // Basic simulation parameter and other misc. channels
    this->time_series_row_buffer_[this->index_map_.time_seconds] =
        static_cast<double>(time_step) * this->parameters.h;
    this->time_series_row_buffer_[this->index_map_.num_convergence_iterations] =
        static_cast<double>(this->solver.convergence_err.size());
    this->time_series_row_buffer_[this->index_map_.convergence_error] =
        this->solver.convergence_err.empty() ? 0. : this->solver.convergence_err.back();
    this->time_series_row_buffer_[this->index_map_.azimuth_angle_degrees] =
        this->CalculateAzimuthAngle() * rad_to_deg;
    this->time_series_row_buffer_[this->index_map_.rotor_speed_rpm] =
        this->CalculateRotorSpeed() / rpm_to_radps;
    this->time_series_row_buffer_[this->index_map_.yaw_position_degrees] =
        this->turbine.yaw_control * rad_to_deg;

    // Tower top and base state data (3 values each)
    for (size_t i = 0; i < 3; ++i) {
        this->time_series_row_buffer_[this->index_map_.tower_top_displacement_start + i] =
            this->turbine.yaw_bearing_node.displacement[i];
        this->time_series_row_buffer_[this->index_map_.tower_top_velocity_start + i] =
            this->turbine.yaw_bearing_node.velocity[i];
        this->time_series_row_buffer_[this->index_map_.tower_top_acceleration_start + i] =
            this->turbine.yaw_bearing_node.acceleration[i];
        this->time_series_row_buffer_[this->index_map_.tower_base_force_start + i] =
            this->turbine.tower_base.loads[i] / 1000.;
        this->time_series_row_buffer_[this->index_map_.tower_base_moment_start + i] =
            this->turbine.tower_base.loads[i + 3] / 1000.;
    }

    // Rotor thrust
    {
        const auto& position = this->turbine.azimuth_node.position;
        const auto q_global_local =
            Eigen::Quaternion<double>(position[3], position[4], position[5], position[6]).inverse();
        const auto shaft_loads =
            Eigen::Matrix<double, 3, 1>(this->turbine.shaft_base_to_azimuth.loads.data());
        const auto shaft_forces = q_global_local._transformVector(shaft_loads);
        this->time_series_row_buffer_[this->index_map_.rotor_thrust_kN] = shaft_forces[0] / 1000.;
    }

    // Blade data
    for (size_t i = 0; i < this->turbine.blades.size(); ++i) {
        const size_t base = this->index_map_.blade_channel_offsets[i];

        const auto& position = this->turbine.blades[i].nodes[0].position;
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
        const auto blade_root_moments = q_global_to_local._transformVector(
            Eigen::Matrix<double, 3, 1>(&this->turbine.blade_pitch[i].loads[3])
        );

        this->time_series_row_buffer_[base + 0] = blade_root_forces[0];
        this->time_series_row_buffer_[base + 1] = blade_root_forces[1];
        this->time_series_row_buffer_[base + 2] = blade_root_forces[2];
        this->time_series_row_buffer_[base + 3] = blade_root_moments[0];
        this->time_series_row_buffer_[base + 4] = blade_root_moments[1];
        this->time_series_row_buffer_[base + 5] = blade_root_moments[2];
        this->time_series_row_buffer_[base + 6] = this->turbine.blade_pitch_control[i] * rad_to_deg;

        const auto& tip = this->turbine.blades[i].nodes.back();
        this->time_series_row_buffer_[base + 7] = tip.velocity[0];
        this->time_series_row_buffer_[base + 8] = tip.velocity[1];
        this->time_series_row_buffer_[base + 9] = tip.velocity[2];
        this->time_series_row_buffer_[base + 10] = tip.velocity[3] / deg_to_rad;
        this->time_series_row_buffer_[base + 11] = tip.velocity[4] / deg_to_rad;
        this->time_series_row_buffer_[base + 12] = tip.velocity[5] / deg_to_rad;
    }

    // Controller data
    if (this->index_map_.has_controller_channels && this->controller) {
        this->time_series_row_buffer_[this->index_map_.generator_torque_kNm] =
            this->controller->io.generator_torque_command / 1000.;
        this->time_series_row_buffer_[this->index_map_.generator_power_kW] =
            this->controller->io.generator_power_actual / 1000.;
    }

    // Hub inflow
    for (size_t i = 0; i < 3; ++i) {
        this->time_series_row_buffer_[this->index_map_.hub_inflow_start + i] = this->hub_inflow[i];
    }

    // Aerodynamic data
    if (this->aerodynamics && !this->index_map_.aero_body_offsets.empty()) {
        const size_t n_bodies = this->index_map_.aero_body_offsets.size();
        for (size_t i = 0; i < n_bodies; ++i) {
            const auto& body = this->aerodynamics->bodies[i];
            const size_t body_base = this->index_map_.aero_body_offsets[i];
            const size_t stride = TimeSeriesIndexMap::kAeroChannelStride;

            for (size_t j = 0; j < body.loads.size(); ++j) {
                const size_t base = body_base + j * stride;

                this->time_series_row_buffer_[base + 0] =
                    Eigen::Matrix<double, 3, 1>(body.v_rel[j].data()).norm();
                this->time_series_row_buffer_[base + 1] = body.alpha[j] / deg_to_rad;
                this->time_series_row_buffer_[base + 2] = body.cn[j];
                this->time_series_row_buffer_[base + 3] = body.ct[j];
                this->time_series_row_buffer_[base + 4] = body.cm[j];
                this->time_series_row_buffer_[base + 5] = body.loads[j][0] / body.delta_s[j];
                this->time_series_row_buffer_[base + 6] = body.loads[j][1] / body.delta_s[j];
                this->time_series_row_buffer_[base + 7] = body.loads[j][2] / body.delta_s[j];
                this->time_series_row_buffer_[base + 8] = body.loads[j][3] / body.delta_s[j];
                this->time_series_row_buffer_[base + 9] = body.loads[j][4] / body.delta_s[j];
                this->time_series_row_buffer_[base + 10] = body.loads[j][5] / body.delta_s[j];
            }
        }
    }

    // write the time-series row at the current time step
    this->outputs->WriteTimeSeriesRowAtTimestep(time_step, this->time_series_row_buffer_);
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

void TurbineInterface::SetHubInflow(const std::array<double, 3>& inflow) {
    this->hub_inflow = inflow;
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
    const double generator_torque =
        this->turbine.torque_control / this->turbine.GetTurbineInput().gear_box_ratio;
    controller->io.horizontal_wind_speed = sqrt(
        (this->hub_inflow[0] * this->hub_inflow[0]) + (this->hub_inflow[1] * this->hub_inflow[1]) +
        (this->hub_inflow[2] * this->hub_inflow[2])
    );
    controller->io.generator_torque_actual = generator_torque;
    controller->io.generator_power_actual =
        generator_speed * generator_torque * this->turbine.GetTurbineInput().generator_efficiency;
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

    this->turbine.torque_control =
        controller->io.generator_torque_command * this->turbine.GetTurbineInput().gear_box_ratio;
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
