#pragma once

#include "interfaces/components/beam.hpp"
#include "interfaces/components/beam_input.hpp"
#include "interfaces/components/outputs_config.hpp"
#include "interfaces/components/solution_input.hpp"
#include "interfaces/host_state.hpp"
#include "interfaces/outputs.hpp"
#include "model/model.hpp"

namespace kynema::interfaces {

/**
 * @brief Interface for blade simulation that manages state, solver, and components
 *
 * This class represents the primary interface for simulating a WT blade, connecting
 * the blade components with the solver and state management.
 */
class BladeInterface {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /**
     * @brief Constructs a BladeInterface from solution and blade inputs
     * @param solution_input Configuration parameters for solver and solution
     * @param blade_input Configuration parameters for blade geometry
     * @param outputs_config Configuration parameters for the outputs
     */
    explicit BladeInterface(
        const components::SolutionInput& solution_input, const components::BeamInput& blade_input,
        const components::OutputsConfig& outputs_config
    );

    /// @brief Returns a reference to the blade model
    [[nodiscard]] components::Beam& Blade();

    /**
     * @brief Steps forward in time
     * @return true if solver converged, false otherwise
     */
    [[nodiscard]] bool Step();

    /// @brief Saves the current state for potential restoration (in correction step)
    void SaveState();

    /// @brief Restores the previously saved state (in correction step)
    void RestoreState();

    /**
     * @brief Sets the displacement for the root node
     * @param u Displacement array (7 components)
     * @throws std::runtime_error if prescribed root motion was not enabled
     */
    void SetRootDisplacement(const std::array<double, 7>& u) const;

private:
    Model model;                    ///< Kynema class for model construction
    components::Beam blade;         ///< Blade model input/output data
    State<DeviceType> state;        ///< Kynema class for storing system state
    Elements<DeviceType> elements;  ///< Kynema class for model elements (beams, masses, springs)
    Constraints<DeviceType> constraints;  ///< Kynema class for constraints tying elements together
    StepParameters parameters;            ///< Kynema class containing solution parameters
    Solver<DeviceType> solver;            ///< Kynema class for solving the dynamic system
    State<DeviceType> state_save;         ///< Kynema class state class for temporarily saving state
    HostState<DeviceType> host_state;     ///< Host local copy of node state data
    std::unique_ptr<Outputs> outputs;     ///< handle to Output for writing to NetCDF
};

}  // namespace kynema::interfaces
