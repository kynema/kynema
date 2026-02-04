#pragma once

#include <functional>
#include <string>

#include "controller_input.hpp"
#include "controller_io.hpp"
#include "vendor/dylib/dylib.hpp"

namespace kynema::interfaces::components {

/// A turbine controller class that works as a wrapper around the shared library containing the
/// controller logic
class Controller {
public:
    /// Pointer to structure mapping swap array -> named fields i.e. ControllerIO
    ControllerIO io;

    /// @brief Constructor for the Controller class
    /// @param shared_lib_path Path to the shared library containing the controller function
    /// @param controller_function_name Name of the controller function in the shared library
    /// @param input_file_path Path to the input file
    /// @param output_file_path Path to the output file
    /// @param initial_yaw_angle Initial yaw angle (rad)
    /// @param yaw_control_enabled Flag to enable yaw control
    Controller(const ControllerInput& input);

    /// Method to call the controller function from the shared library
    void CallController();

    /// @brief Get the commanded pitch angle (rad)
    [[nodiscard]] double PitchAngleCommand() const { return this->pitch_angle_command_; }

    /// @brief Get the commanded torque (Nm)
    [[nodiscard]] double TorqueCommand() const { return this->torque_command_; }

    /// @brief Get the commanded yaw angle (rad) integrated from yaw rate command
    [[nodiscard]] double YawAngleCommand() const { return this->yaw_angle_command_; }

private:
    bool pitch_control_enabled_;   //< Flag to enable pitch control
    bool torque_control_enabled_;  //< Flag to enable torque control
    bool yaw_control_enabled_;     //< Flag to enable yaw control

    double gearbox_ratio_;  // Gearbox ratio

    double torque_command_;       //< Commanded torque (Nm)
    double pitch_angle_command_;  //< Commanded pitch angle (rad)
    double yaw_angle_command_;    //< Commanded yaw angle (rad) integrated from yaw rate command

    std::string input_file_path_;           //< Path to the input file
    std::string output_file_path_;          //< Path to the output file
    std::string shared_lib_path_;           //< Path to shared library
    std::string controller_function_name_;  //< Name of the controller function in the shared library

    util::dylib lib_;  //< Handle to the shared library
    std::function<void(float*, int*, const char* const, char* const, char* const)>
        controller_function_;  //< Function pointer to the controller function
};

}  // namespace kynema::interfaces::components
