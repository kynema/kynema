#pragma once

#include <string>
#include <string_view>

#include "controller_input.hpp"

namespace kynema::interfaces::components {

class ControllerBuilder {
public:
    ControllerInput& Input() { return input; }

    // Enable or disable the controller
    ControllerBuilder& EnableController(bool enable = true) {
        input.controller_enabled = enable;
        return *this;
    }

    // Enable or disable pitch control
    ControllerBuilder& EnablePitchControl(bool enable = true) {
        input.pitch_control_enabled = enable;
        return *this;
    }

    // Enable or disable torque control
    ControllerBuilder& EnableTorqueControl(bool enable = true) {
        input.torque_control_enabled = enable;
        return *this;
    }

    // Set time step (seconds)
    ControllerBuilder& SetTimeStep(double time_step) {
        input.time_step = time_step;
        return *this;
    }

    // Set number of blades
    ControllerBuilder& SetNumberOfBlades(uint n_blades) {
        input.n_blades = n_blades;
        return *this;
    }

    // Enable or disable yaw control
    ControllerBuilder& EnableYawControl(bool enable = true) {
        input.yaw_control_enabled = enable;
        return *this;
    }

    // Set initial yaw angle (radians)
    ControllerBuilder& SetYawAngle(double angle) {
        input.yaw_angle = angle;
        return *this;
    }

    // Set initial rotor speed (rad/s)
    ControllerBuilder& SetRotorSpeed(double speed) {
        input.rotor_speed = speed;
        return *this;
    }

    // Set initial power (W)
    ControllerBuilder& SetPower(double power) {
        input.power = power;
        return *this;
    }

    /// @brief Set flag to enable reading of checkpoint file for restart
    ControllerBuilder& EnableReadCheckpoint(bool enable = true) {
        input.read_checkpoint = enable;
        return *this;
    }

    /// @brief Set library path for the controller
    ControllerBuilder& SetLibraryPath(std::string_view lib_path) {
        input.shared_lib_path = std::string(lib_path);
        return *this;
    }

    /// @brief Set function name for the controller
    ControllerBuilder& SetFunctionName(std::string_view func_name) {
        input.function_name = std::string(func_name);
        return *this;
    }

    /// @brief Set input file path for the controller
    ControllerBuilder& SetInputFilePath(std::string_view inp_file_path) {
        input.input_file_path = std::string(inp_file_path);
        return *this;
    }

    /// @brief Set output file path for the controller
    ControllerBuilder& SetOutputFilePath(std::string_view sim_name) {
        input.output_file_path = std::string(sim_name);
        return *this;
    }

private:
    ControllerInput input;
};

}  // namespace kynema::interfaces::components
