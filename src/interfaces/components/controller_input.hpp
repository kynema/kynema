#pragma once

#include <sys/types.h>
#include <string>
#include <cstdint>

namespace kynema::interfaces::components {

enum class PitchActuatorType : uint8_t {
    POSITION = 0,
    VELOCITY = 1,
};

enum class PitchControlType : uint8_t {
    COLLECTIVE = 0,
    INDIVIDUAL = 1,
};

/**
 * @brief Configuration parameters for a DISCON-style turbine controller
 *
 * This struct encapsulates all the necessary configuration parameters needed to
 * initialize and configure a DISCON-style wind turbine controller. DISCON is a
 * standardized interface for wind turbine control algorithms that allows for
 * dynamic loading of controller implementations at runtime.
 */
struct ControllerInput {
    bool controller_enabled{false};      ///< Flag to enable controller
    bool pitch_control_enabled{false};   ///< Flag to enable pitch control
    bool torque_control_enabled{false};  ///< Flag to enable torque control
    bool yaw_control_enabled{false};     ///< Flag to enable yaw control
    PitchActuatorType pitch_actuator_type{PitchActuatorType::POSITION
    };  ///< Pitch actuator type (default: 0 for position)
    PitchControlType pitch_control_type{PitchControlType::COLLECTIVE
    };  ///< Pitch control type (default: 0 for collective)
    double time_step{0.0};         ///< Time step for controller simulation (s)
    uint n_blades{0};              ///< Number of blades
    double pitch_angle{0.0};       ///< Initial pitch angle (radians)
    double yaw_angle{0.0};         ///< Initial yaw angle (radians)
    double rotor_speed{0.0};       ///< Initial rotor speed (rad/s)
    double power{0.0};             ///< Initial power (W)
    std::string shared_lib_path;   ///< Path to controller shared library
    std::string function_name;     ///< Controller function name (default: "DISCON")
    std::string input_file_path;   ///< Path to controller input file
    std::string output_file_path;  ///< Simulation name for controller
};

}  // namespace kynema::interfaces::components
