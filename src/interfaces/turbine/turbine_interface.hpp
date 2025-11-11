#pragma once

#include "interfaces/components/aerodynamics.hpp"
#include "interfaces/components/aerodynamics_input.hpp"
#include "interfaces/components/controller_input.hpp"
#include "interfaces/components/outputs_config.hpp"
#include "interfaces/components/turbine.hpp"
#include "interfaces/host_constraints.hpp"
#include "interfaces/host_state.hpp"
#include "interfaces/outputs.hpp"
#include "model/model.hpp"
#include "step/step_parameters.hpp"
#include "utilities/controllers/turbine_controller.hpp"

namespace kynema::interfaces::components {
struct SolutionInput;
struct TurbineInput;
struct OutputsConfig;
}  // namespace kynema::interfaces::components

namespace kynema::interfaces {

/**
 * @brief Interface for blade simulation that manages state, solver, and components
 *
 * This class represents the primary interface for simulating a WT blade, connecting
 * the blade components with the solver and state management.
 */
class TurbineInterface {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /**
     * @brief Constructs a TurbineInterface from solution and blade inputs
     * @param solution_input Configuration parameters for solver and solution
     * @param turbine_input Configuration parameters for the turbine geometry
     * @param aerodynamics_input Configuration parameters for the aerodynamic loads
     * @param controller_input Configuration parameters for the controller
     * @param outputs_config Configuration parameters for the outputs
     */
    explicit TurbineInterface(
        const components::SolutionInput& solution_input,
        const components::TurbineInput& turbine_input,
        const components::AerodynamicsInput& aerodynamics_input = {},
        const components::ControllerInput& controller_input = {},
        const components::OutputsConfig& outputs_config = {}
    );

    /// @brief Returns a reference to the turbine model
    [[nodiscard]] components::Turbine& Turbine() { return this->turbine; }

    /// @brief Returns a reference to the aerodynamics model
    [[nodiscard]] components::Aerodynamics& Aerodynamics() {
        if (!aerodynamics) {
            throw std::runtime_error("Aerodynamics component not initialized in TurbineInterface.");
        }
        return *aerodynamics;
    }

    /**
     * @brief Updates the aerodynamic loads to be applied to the structure based on a provided
     * function
     *
     * @param fluid_density The density of the air (assumed constant)
     * @param inflow_function A function that takes in a position and returns the velocity
     */
    void UpdateAerodynamicLoads(
        double fluid_density,
        const std::function<std::array<double, 3>(const std::array<double, 3>&)>& inflow_function
    );

    std::array<double, 3> GetHubNodePosition() const;

    /**
     * @brief Update controller inputs from current system state
     */
    void ApplyController(double t);

    /**
     * @brief Steps forward in time
     *
     * @return true if solver converged, false otherwise
     * @note This function updates the host state with current node loads,
     *       solves the dynamic system, and updates the node motion with the new state.
     *       If the solver does not converge, the motion is not updated.
     */
    [[nodiscard]] bool Step();

    /// @brief Saves the current state for potential restoration (in correction step)
    void SaveState();

    /// @brief Restores the previously saved state (in correction step)
    void RestoreState();

    /// @brief Return a reference of the model owned by this interface
    Model& GetModel() { return model; }

    /// @brief Return a reference to this interface's host state
    HostState<DeviceType>& GetHostState() { return host_state; }

    /**
     * @brief Calculates and normalizes azimuth angle from constraint output
     * @return Azimuth angle in radians, normalized to [0, 2π)
     */
    [[nodiscard]] double CalculateAzimuthAngle() const;

    /**
     * @brief Calculates rotor speed from constraint output
     * @return Rotor speed in rad/s
     */
    [[nodiscard]] double CalculateRotorSpeed() const;

    void WriteOutput();

    void OpenOutputFile();

    void CloseOutputFile();

private:
    Model model;                    ///< Kynema class for model construction
    components::Turbine turbine;    ///< Turbine model input/output data
    State<DeviceType> state;        ///< Kynema class for storing system state
    Elements<DeviceType> elements;  ///< Kynema class for model elements (beams, masses, springs)
    Constraints<DeviceType> constraints;  ///< Kynema class for constraints tying elements together
    StepParameters parameters;            ///< Kynema class containing solution parameters
    Solver<DeviceType> solver;            ///< Kynema class for solving the dynamic system
    State<DeviceType> state_save;         ///< Kynema class state class for temporarily saving state
    HostState<DeviceType> host_state;     ///< Host local copy of node state data
    HostConstraints<DeviceType> host_constraints;         ///< Host local copy of constraint data
    std::unique_ptr<Outputs> outputs;                     ///< handle to Output for writing to NetCDF
    std::unique_ptr<util::TurbineController> controller;  ///< DISCON-style controller
    std::unique_ptr<components::Aerodynamics> aerodynamics;  ///< Aerodynamics component
    std::array<double, 3> hub_inflow{0., 0., 0.};            ///< Inflow velocity at the hub node

    /**
     * @brief Write rotor time-series data based on constraint outputs
     *
     * This method extracts rotor azimuth angle and speed from the constraint system
     * and writes them to the time-series output file. Data is read from the shaft
     * base to azimuth constraint, which contains:
     * - Index 0: Azimuth angle (radians)
     * - Index 1: Rotor speed (rad/s)
     */
    void WriteTimeSeriesData() const;

    /**
     * @brief Initialize controller with turbine parameters and connect to constraints
     * @param turbine_input Configuration parameters for turbine geometry and initial conditions
     */
    void InitializeController(
        const components::TurbineInput& turbine_input,
        const components::SolutionInput& solution_input
    );
};

}  // namespace kynema::interfaces
