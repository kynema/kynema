#include "turbine_interface.hpp"

#include <filesystem>
#include <numbers>

#include "interfaces/components/solution_input.hpp"
#include "state/clone_state.hpp"
#include "state/copy_state_data.hpp"
#include "step/step.hpp"

namespace kynema::interfaces {

TurbineInterface::TurbineInterface(
    const components::SolutionInput& solution_input, const components::TurbineInput& turbine_input,
    const components::AerodynamicsInput& aerodynamics_input,
    const components::ControllerInput& controller_input
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
      host_state(state) {
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

    // Initialize NetCDF writer and write mesh connectivity if output path is specified
    if (!solution_input.output_file_path.empty()) {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(solution_input.output_file_path);

        // Write mesh connectivity to YAML file
        model.ExportMeshConnectivityToYAML(
            solution_input.output_file_path + "/mesh_connectivity.yaml"
        );

        // Initialize outputs with both node state and time-series files
        this->outputs = std::make_unique<Outputs>(
            solution_input.output_file_path + "/turbine_interface.nc",
            solution_input.output_file_path + "/turbine_time_series.nc", this->state.num_system_nodes
        );

        // Write initial state
        this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

        // Write initial time-series data (test values)
        this->outputs->WriteRotorTimeSeriesAtTimestep(
            this->state.time_step, turbine_input.azimuth_angle, turbine_input.rotor_speed
        );
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
        aerodynamics->CalculateMotion(host_state);

        aerodynamics->SetInflowFromFunction(inflow_function);

        aerodynamics->CalculateAerodynamicLoads(fluid_density);

        aerodynamics->CalculateNodalLoads();
    }
}

bool TurbineInterface::Step(double t) {
    // Apply controller if available
    if (controller) {
        controller->io.time = t;
        ApplyController();
    }

    // Update the host state with current node loads
    Kokkos::deep_copy(this->host_state.f, 0.);
    this->turbine.SetLoads(this->host_state);
    if (this->aerodynamics) {
        this->aerodynamics->AddNodalLoadsToState(this->host_state);
    }
    Kokkos::deep_copy(this->state.f, this->host_state.f);

    // Solve for state at end of step
    auto converged =
        kynema::Step(this->parameters, this->solver, this->elements, this->state, this->constraints);

    // If not converged, return false
    if (!converged) {
        return false;
    }

    // Update the host state with current node motion
    this->host_state.CopyFromState(this->state);

    // Update the turbine node motion based on the host state
    this->turbine.GetMotion(this->host_state);

    // Write outputs and increment timestep counter
    if (this->outputs) {
        // Write node state outputs
        this->outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);

        // Calculate rotor azimuth and speed -> write rotor time-series data
        this->WriteRotorTimeSeriesData();
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

void TurbineInterface::WriteRotorTimeSeriesData() {
    if (!this->outputs) {
        return;
    }

    // Write the calculated values to time-series output
    this->outputs->WriteRotorTimeSeriesAtTimestep(
        this->state.time_step,
        this->CalculateAzimuthAngle(),  // azimuth angle (radians)
        this->CalculateRotorSpeed()     // rotor speed (rad/s)
    );
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
}

void TurbineInterface::ApplyController() {
    if (!controller) {
        return;
    }

    // Update controller inputs from current system state
    UpdateControllerInputs();

    // Call the controller
    controller->CallController();

    // Apply controller outputs to the turbine?
}

void TurbineInterface::UpdateControllerInputs() {
    // Update time and azimuth
    // TODO How to get simulation time here?
    controller->io.time = static_cast<double>(this->state.time_step);
    controller->io.azimuth_angle = CalculateAzimuthAngle();

    // Update rotor and generator speeds
    const double rotor_speed = CalculateRotorSpeed();
    controller->io.rotor_speed_actual = rotor_speed;
    controller->io.generator_speed_actual =
        rotor_speed * this->turbine.GetTurbineInput().gear_box_ratio;

    // Update generator power and torque
    const double generator_speed = controller->io.generator_speed_actual;
    const double generator_torque = controller->io.generator_torque_actual;
    controller->io.generator_power_actual = generator_speed * generator_torque;

    // Update wind speed (assuming it's constant for now)
    controller->io.horizontal_wind_speed = this->turbine.GetTurbineInput().hub_wind_speed;

    // Set status for subsequent calls
    controller->io.status = 1;
}
}  // namespace kynema::interfaces
