#pragma once

#include <array>

namespace kynema::util {

static constexpr int kSwapArraySize = 128;

/// Structure defining the swap layout (See Appendix A of Bladed User's Guide)
/// Ref: https://openfast.readthedocs.io/en/main/source/user/servodyn/ExtendedBladedInterface.html
struct ControllerIO {
    int status = 0;    // Input: Status flag [0=first call, 1=intermediate, -1=final] (-)
    double time = 0.;  // Input: Current time (sec)
    double dt = 0.;    // Input: Communication interval (sec)
    double pitch_blade1_actual = 0.;      // Input: Blade 1 pitch angle (rad)
    double pitch_actuator_type_req = 0.;  // Input: Pitch actuator type [0=position, 1=rate] (-)
    double generator_power_actual = 0.;   // Input: Measured electrical power output (W)
    double generator_speed_actual = 0.;   // Input: Measured generator speed (rad/s)
    double rotor_speed_actual = 0.;       // Input: Measured rotor speed (rad/s)
    double generator_torque_actual = 0.;  // Input: Measured generator torque (Nm)
    double yaw_error_actual = 0.;         // Input: Measured yaw error (rad)
    double horizontal_wind_speed = 0.;    // Input: Horizontal hub-height wind speed (m/s)
    double pitch_control_type = 0.;       // Input: Pitch control type [0=collective, 1=individual]
    double out_of_plane_root_bending_moment_blade1 =
        0.;  // Input: Out-of-plane root bending moment blade 1 (Nm)
    double out_of_plane_root_bending_moment_blade2 =
        0.;  // Input: Out-of-plane root bending moment blade 2 (Nm)
    double out_of_plane_root_bending_moment_blade3 =
        0.;                                   // Input: Out-of-plane root bending moment blade 3 (Nm)
    double pitch_blade2_actual = 0.;          // Input: Blade 2 pitch
    double pitch_blade3_actual = 0.;          // Input: Blade 3 pitch
    int generator_contactor_status = 0;       // Output: Generator contactor status
    int shaft_brake_status = 0;               // Input: Shaft brake status [0=off, 1=on]
    double yaw_angle_actual = 0;              // Input: Nacelle yaw angle from North (rad)
    double yaw_actuator_torque_command = 0.;  // Output: Demanded yaw actuator torque
    double pitch_blade1_command = 0.;         // Output: pitch command (rad)
    double pitch_blade2_command = 0.;         // Output: pitch command (rad)
    double pitch_blade3_command = 0.;         // Output: pitch command (rad)
    double pitch_collective_command = 0.;     // Output: pitch collective command (rad)
    double pitch_rate_command = 0.;           // Output: Demanded pitch rate (Collective pitch)
    double generator_torque_command = 0.;     // Output: Demanded generator torque
    double nacelle_yaw_rate_command = 0.;     // Output: Demanded nacelle yaw rate
    size_t message_array_size = 0U;           // Input: avcMSG array size
    size_t infile_array_size = 0U;            // Input: avcINFILE array size
    size_t outname_array_size = 0U;           // Input: avcOUTNAME array size
    double tower_top_fore_aft_accel = 0.;     // Input: Tower top fore-aft  acceleration (m/s^2)
    int pitch_override = 0;                   // Input: Pitch override: 0=yes
    int torque_override = 0;                  // Input: Torque override: 0=yes
    double azimuth_angle = 0.;                // Input: Azimuth angle (rad)
    size_t n_blades = 0U;                     // Input: Number of blades
    size_t n_log_variables = 0U;              // Input: Number of variables returned for logging
    double generator_startup_resistance = 0.;  // Input: Generator start-up resistance
    int loads_request = 0;                     // Input: Request for loads: 0=none
    int variable_slip_status = 0;              // Input: Variable slip current status
    int variable_slip_demand = 0;              // Input: Variable slip current demand
    double nacelle_nodding_accel = 0.;         // Input: Nacelle nodding acceleration (rad/s^2)

    void CopyToSwapArray(std::array<float, kSwapArraySize>& swap_array) const {
        swap_array[0] = static_cast<float>(status);
        swap_array[1] = static_cast<float>(time);
        swap_array[2] = static_cast<float>(dt);
        swap_array[3] = static_cast<float>(pitch_blade1_actual);
        swap_array[9] = static_cast<float>(pitch_actuator_type_req);
        swap_array[14] = static_cast<float>(generator_power_actual);
        swap_array[19] = static_cast<float>(generator_speed_actual);
        swap_array[20] = static_cast<float>(rotor_speed_actual);
        swap_array[22] = static_cast<float>(generator_torque_actual);
        swap_array[23] = static_cast<float>(yaw_error_actual);
        swap_array[26] = static_cast<float>(horizontal_wind_speed);
        swap_array[27] = static_cast<float>(pitch_control_type);
        swap_array[29] = static_cast<float>(out_of_plane_root_bending_moment_blade1);
        swap_array[30] = static_cast<float>(out_of_plane_root_bending_moment_blade2);
        swap_array[31] = static_cast<float>(out_of_plane_root_bending_moment_blade3);
        swap_array[32] = static_cast<float>(pitch_blade2_actual);
        swap_array[33] = static_cast<float>(pitch_blade3_actual);
        swap_array[34] = static_cast<float>(generator_contactor_status);
        swap_array[35] = static_cast<float>(shaft_brake_status);
        swap_array[36] = static_cast<float>(yaw_angle_actual);
        swap_array[40] = static_cast<float>(yaw_actuator_torque_command);
        swap_array[41] = static_cast<float>(pitch_blade1_command);
        swap_array[42] = static_cast<float>(pitch_blade2_command);
        swap_array[43] = static_cast<float>(pitch_blade3_command);
        swap_array[44] = static_cast<float>(pitch_collective_command);
        swap_array[45] = static_cast<float>(pitch_rate_command);
        swap_array[46] = static_cast<float>(generator_torque_command);
        swap_array[47] = static_cast<float>(nacelle_yaw_rate_command);
        swap_array[48] = static_cast<float>(message_array_size);
        swap_array[49] = static_cast<float>(infile_array_size);
        swap_array[50] = static_cast<float>(outname_array_size);
        swap_array[52] = static_cast<float>(tower_top_fore_aft_accel);
        swap_array[54] = static_cast<float>(pitch_override);
        swap_array[55] = static_cast<float>(torque_override);
        swap_array[59] = static_cast<float>(azimuth_angle);
        swap_array[60] = static_cast<float>(n_blades);
        swap_array[64] = static_cast<float>(n_log_variables);
        swap_array[71] = static_cast<float>(generator_startup_resistance);
        swap_array[78] = static_cast<float>(loads_request);
        swap_array[79] = static_cast<float>(variable_slip_status);
        swap_array[80] = static_cast<float>(variable_slip_demand);
        swap_array[82] = static_cast<float>(nacelle_nodding_accel);
    }

    void CopyFromSwapArray(const std::array<float, kSwapArraySize>& swap_array) {
        status = static_cast<int>(swap_array[0]);
        time = static_cast<double>(swap_array[1]);
        dt = static_cast<double>(swap_array[2]);
        pitch_blade1_actual = static_cast<double>(swap_array[3]);
        pitch_actuator_type_req = static_cast<double>(swap_array[9]);
        generator_power_actual = static_cast<double>(swap_array[14]);
        generator_speed_actual = static_cast<double>(swap_array[19]);
        rotor_speed_actual = static_cast<double>(swap_array[20]);
        generator_torque_actual = static_cast<double>(swap_array[22]);
        yaw_error_actual = static_cast<double>(swap_array[23]);
        horizontal_wind_speed = static_cast<double>(swap_array[26]);
        pitch_control_type = static_cast<int>(swap_array[27]);
        out_of_plane_root_bending_moment_blade1 = static_cast<double>(swap_array[29]);
        out_of_plane_root_bending_moment_blade2 = static_cast<double>(swap_array[30]);
        out_of_plane_root_bending_moment_blade3 = static_cast<double>(swap_array[31]);
        pitch_blade2_actual = static_cast<double>(swap_array[32]);
        pitch_blade3_actual = static_cast<double>(swap_array[33]);
        generator_contactor_status = static_cast<int>(swap_array[34]);
        shaft_brake_status = static_cast<int>(swap_array[35]);
        yaw_angle_actual = static_cast<double>(swap_array[36]);
        yaw_actuator_torque_command = static_cast<double>(swap_array[40]);
        pitch_blade1_command = static_cast<double>(swap_array[41]);
        pitch_blade2_command = static_cast<double>(swap_array[42]);
        pitch_blade3_command = static_cast<double>(swap_array[43]);
        pitch_collective_command = static_cast<double>(swap_array[44]);
        pitch_rate_command = static_cast<double>(swap_array[45]);
        generator_torque_command = static_cast<double>(swap_array[46]);
        nacelle_yaw_rate_command = static_cast<double>(swap_array[47]);
        message_array_size = static_cast<size_t>(swap_array[48]);
        infile_array_size = static_cast<size_t>(swap_array[49]);
        outname_array_size = static_cast<size_t>(swap_array[50]);
        tower_top_fore_aft_accel = static_cast<double>(swap_array[52]);
        pitch_override = static_cast<int>(swap_array[54]);
        torque_override = static_cast<int>(swap_array[55]);
        azimuth_angle = static_cast<double>(swap_array[59]);
        n_blades = static_cast<size_t>(swap_array[60]);
        n_log_variables = static_cast<size_t>(swap_array[64]);
        generator_startup_resistance = static_cast<double>(swap_array[71]);
        loads_request = static_cast<int>(swap_array[78]);
        variable_slip_status = static_cast<int>(swap_array[79]);
        variable_slip_demand = static_cast<int>(swap_array[80]);
        nacelle_nodding_accel = static_cast<double>(swap_array[82]);
    }
};

}  // namespace kynema::util
