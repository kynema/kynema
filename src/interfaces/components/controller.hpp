#pragma once

#include <array>
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

    /// @brief Get the collective commanded pitch angle (rad)
    [[nodiscard]] double PitchAngleCommand() const { return this->pitch_angle_command_; }

    /// @brief Get the commanded individual pitch angle (rad)
    [[nodiscard]] std::array<double, 3> PitchAngleCommandIndividual() const {
        return {
            this->io.pitch_blade1_command, this->io.pitch_blade2_command,
            this->io.pitch_blade3_command
        };
    }

    /// @brief Get the commanded generator torque (Nm)
    [[nodiscard]] double GeneratorTorqueCommand() const { return this->generator_torque_command_; }

    /// @brief Get the commanded yaw angle (rad) integrated from yaw rate command
    [[nodiscard]] double YawAngleCommand() const { return this->yaw_angle_command_; }

    /// @brief Set the controller status to initialization
    void SetStatusInit() { this->io.status = 0; }

    /// @brief Set the controller status to operating
    void SetStatusOperating() { this->io.status = 1; }

    /// @brief Set the controller status to write checkpoint file (only valid for ROSCO controller)
    void SetStatusWriteCheckpoint() { this->io.status = -8; }

    /// @brief Set the controller status to read checkpoint file (only valid for ROSCO controller)
    void SetStatusReadCheckpoint() { this->io.status = -9; }

    /// @brief Set the simulation time (s)
    void SetSimulationTime(double simulation_time) { this->io.time = simulation_time; }

    /// @brief Set the simulation time step (s)
    void SetSimulationTimeStep(double simulation_time_step) { this->io.dt = simulation_time_step; }

    /// @brief Set the measured yaw angle (rad)
    void SetYawAngle(double yaw_angle) { this->io.yaw_angle_actual = yaw_angle; }

    /// @brief Set the measured rotor speed (rad/s)
    void SetRotorSpeed(double rotor_speed) { this->io.rotor_speed_actual = rotor_speed; }

    /// @brief Set the measured generator speed (rad/s)
    void SetGeneratorSpeed(double generator_speed) {
        this->io.generator_speed_actual = generator_speed;
    }

    /// @brief Set the measured generator torque (Nm)
    void SetGeneratorTorque(double generator_torque) {
        this->io.generator_torque_actual = generator_torque;
    }

    /// @brief Set measured generator power (W)
    void SetGeneratorPower(double generator_power) {
        this->io.generator_power_actual = generator_power;
    }

    /// @brief Set measured wind speed at hub (m/s)
    void SetWindSpeed(double wind_speed) { this->io.horizontal_wind_speed = wind_speed; }

    /// @brief Set the rotor azimuth angle (rad)
    void SetRotorAzimuth(double rotor_azimuth) { this->io.azimuth_angle = rotor_azimuth; }

    /// @brief Set measured blade pitch angle (rad)
    void SetBladePitch(double blade_pitch_collective) {
        this->io.pitch_blade1_actual = blade_pitch_collective;
        this->io.pitch_blade2_actual = blade_pitch_collective;
        this->io.pitch_blade3_actual = blade_pitch_collective;
    }

    /// @brief Set measured blade pitch angle (rad)
    void SetBladePitch(std::array<double, 3> blade_pitch) {
        this->io.pitch_blade1_actual = blade_pitch[0];
        this->io.pitch_blade2_actual = blade_pitch[1];
        this->io.pitch_blade3_actual = blade_pitch[2];
    }

    /// @brief Set measured out-of-plane root bending moment for each blade (Nm)
    void SetOutOfPlaneRootBendingMoment(std::array<double, 3> root_moment) {
        this->io.out_of_plane_root_bending_moment_blade1 = root_moment[0];
        this->io.out_of_plane_root_bending_moment_blade2 = root_moment[1];
        this->io.out_of_plane_root_bending_moment_blade3 = root_moment[2];
    }

private:
    /// Pointer to structure mapping swap array -> named fields i.e. ControllerIO
    ControllerIO io;

    bool pitch_control_enabled_;   //< Flag to enable pitch control
    bool torque_control_enabled_;  //< Flag to enable torque control
    bool yaw_control_enabled_;     //< Flag to enable yaw control

    double generator_torque_command_;  //< Commanded torque (Nm)
    double pitch_angle_command_;       //< Commanded pitch angle (rad)
    double yaw_angle_command_;         //< Commanded yaw angle (rad) integrated from yaw rate command

    std::string input_file_path_;           //< Path to the input file
    std::string output_file_path_;          //< Path to the output file
    std::string shared_lib_path_;           //< Path to shared library
    std::string controller_function_name_;  //< Name of the controller function in the shared library

    util::dylib lib_;  //< Handle to the shared library
    std::function<void(float*, int*, const char* const, const char* const, char* const)>
        controller_function_;  //< Function pointer to the controller function
};

}  // namespace kynema::interfaces::components
