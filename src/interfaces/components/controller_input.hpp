#pragma once

#include <string>

namespace kynema::interfaces::components {

/**
 * @brief Configuration parameters for a DISCON-style turbine controller
 *
 * This struct encapsulates all the necessary configuration parameters needed to
 * initialize and configure a DISCON-style wind turbine controller. DISCON is a
 * standardized interface for wind turbine control algorithms that allows for
 * dynamic loading of controller implementations at runtime.
 */
struct ControllerInput {
    std::string shared_lib_path;      ///< Path to controller shared library
    std::string function_name;        ///< Controller function name (default: "DISCON")
    std::string input_file_path;      ///< Path to controller input file
    std::string simulation_name;      ///< Simulation name for controller
    bool yaw_control_enabled{false};  ///< Flag to enable yaw control

    /// @brief Default constructor - creates empty controller input
    ControllerInput() = default;

    /**
     * @brief Constructor with all parameters
     * @param lib_path Path to the shared library containing the controller implementation
     * @param func_name Name of the controller function to call (defaults to "DISCON")
     * @param inp_file_path Optional path to controller-specific configuration file
     * @param sim_name Optional identifier for the simulation run
     */
    explicit ControllerInput(
        std::string lib_path, std::string func_name = "DISCON", std::string inp_file_path = "",
        std::string sim_name = ""
    )
        : shared_lib_path(std::move(lib_path)),
          function_name(std::move(func_name)),
          input_file_path(std::move(inp_file_path)),
          simulation_name(std::move(sim_name)) {}

    /// @brief Check if controller is enabled (i.e. has library path)
    /// @return If the shared library path is set
    [[nodiscard]] bool IsEnabled() const { return !shared_lib_path.empty(); }
};

}  // namespace kynema::interfaces::components
