#pragma once

#include <functional>
#include <string>

#include "utilities/controllers/controller_io.hpp"
#include "vendor/dylib/dylib.hpp"

namespace kynema::util {

/// A turbine controller class that works as a wrapper around the shared library containing the
/// controller logic
class TurbineController {
public:
    /// Pointer to structure mapping swap array -> named fields i.e. ControllerIO
    ControllerIO io;

    /// @brief Constructor for the TurbineController class
    /// @param shared_lib_path Path to the shared library containing the controller function
    /// @param controller_function_name Name of the controller function in the shared library
    /// @param input_file_path Path to the input file
    /// @param output_file_path Path to the output file
    /// @param yaw_control_enabled Flag to enable yaw control
    TurbineController(
        std::string shared_lib_path, std::string controller_function_name,
        std::string input_file_path, std::string output_file_path, bool yaw_control_enabled = false
    );

    /// Method to call the controller function from the shared library
    void CallController();

    /// @brief Get the commanded yaw angle (rad) integrated from yaw rate command
    [[nodiscard]] double YawAngleCommand() const { return this->yaw_angle_command_; }

private:
    std::string input_file_path_;           //< Path to the input file
    std::string output_file_path_;          //< Path to the output file
    std::string shared_lib_path_;           //< Path to shared library
    std::string controller_function_name_;  //< Name of the controller function in the shared library

    util::dylib lib_;  //< Handle to the shared library
    std::function<void(float*, int*, const char* const, char* const, char* const)>
        controller_function_;  //< Function pointer to the controller function

    //< Commanded yaw angle (rad) integrated from yaw rate command
    double yaw_angle_command_{0.0};

    //< Flag to enable yaw control
    bool yaw_control_enabled_{false};
};

}  // namespace kynema::util
