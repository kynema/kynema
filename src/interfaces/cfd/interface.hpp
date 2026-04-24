#pragma once

#include <filesystem>
#include <memory>

#include "constraints/constraints.hpp"
#include "elements/elements.hpp"
#include "interfaces/cfd/interface_input.hpp"
#include "interfaces/cfd/turbine.hpp"
#include "interfaces/host_state.hpp"
#include "interfaces/outputs.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"
#include "state/state.hpp"
#include "step/step_parameters.hpp"

namespace kynema_fmb::interfaces::cfd {

/**
 * @brief The main interface for controlling the CFD problem
 */
class Interface {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /**
     * @brief Create an Interface from an InterfaceInput configuration object
     *
     * @param input The input configuration object
     */
    explicit Interface(const InterfaceInput& input);

    /// @brief Step forward in time one time step
    [[nodiscard]] bool Step();

    /// @brief Save state for correction step
    void SaveState();

    /// @brief Restore state for correction step
    void RestoreState();

    /**
     * @brief Write restart file
     *
     * @param filename The name of the restart file
     */
    void WriteRestart(const std::filesystem::path& filename) const;

    /**
     * @brief Read restart file
     *
     * @param filename The name of the restart file
     */
    void ReadRestart(const std::filesystem::path& filename);

    /// @brief Write current state to output file if configured
    void WriteOutputs() const;

    /// @brief  Kynema class used for model construction
    Model model;

    /// @brief Turbine model input/output data
    Turbine turbine;

    /// @brief  Kynema class for storing system state
    State<DeviceType> state;

    /// @brief  Kynema class for model elements (beams, masses, springs)
    Elements<DeviceType> elements;

    /// @brief  Kynema class for constraints tying elements together
    Constraints<DeviceType> constraints;

    /// @brief  Kynema class containing solution parameters
    StepParameters parameters;

    /// @brief  Kynema class for solving the dynamic system
    Solver<DeviceType> solver;

    /// @brief  Kynema class state class for temporarily saving state
    State<DeviceType> state_save;

    /// @brief Host local copy of State
    kynema_fmb::interfaces::HostState<DeviceType> host_state;

    /// @brief Current timestep index
    size_t current_timestep_{0};

    /// @brief Optional NetCDF output writer
    std::unique_ptr<interfaces::Outputs> outputs_;
};

}  // namespace kynema_fmb::interfaces::cfd
