#include "discon.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <ranges>

#include "controller_io.hpp"

namespace kynema::util {

extern "C" {

// Define some constants
// ------------------------------------------------------------------------------------------------
/// Transitional generator speed (HSS side) between regions 1 and 1 1/2, rad/s
static constexpr double kVS_CtInSp{70.16224};
/// Communication interval for torque controller, sec
static constexpr double kVS_DT{0.000125};
/// Maximum torque rate (in absolute value) in torque controller, N-m/s
static constexpr double kVS_MaxRat{15000.};
/// Maximum generator torque in Region 3 (HSS side), N-m
static constexpr double kVS_MaxTq{47402.91};
/// Generator torque constant in Region 2 (HSS side), N-m/(rad/s)^2
static constexpr double kVS_Rgn2K{2.332287};
/// Transitional generator speed (HSS side) between regions 1 1/2 and 2, rad/s
static constexpr double kVS_Rgn2Sp{91.21091};
/// Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the
/// generator speed, rad
static constexpr double kVS_Rgn3MP{0.01745329};
/// Rated generator speed (HSS side), rad/s
static constexpr double kVS_RtGnSp{121.6805};
/// Rated generator generator power in Region 3, Watts
static constexpr double kVS_RtPwr{5296610.0};
/// Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s -- chosen to
/// be 1/4 the blade edgewise natural frequency ( 1/4 of approx. 1 Hz = 0.25 Hz = 1.570796 rad/s)
static constexpr double kCornerFreq{1.570796};
/// A value slightly greater than unity in single precision
static constexpr double kOnePlusEps{1.0 + 1.19e-07};
/// Communication interval for the pitch controller, sec
static constexpr double kPC_DT{0.000125};
/// Integral gain for pitch controller at rated pitch (zero), (-)
static constexpr double kPC_KI{0.008068634};
/// Pitch angle where the the derivative of the aerodynamic power w.r.t. pitch has increased by a
/// factor of two relative to the derivative at rated pitch (zero), rad
static constexpr double kPC_KK{0.1099965};
/// Proportional gain for pitch controller at rated pitch (zero), sec
static constexpr double kPC_KP{0.01882681};
/// Maximum pitch setting in pitch controller, rad
static constexpr double kPC_MaxPit{1.570796};
/// Maximum pitch rate (in absolute value) in pitch controller, rad/s
static constexpr double kPC_MaxRat{0.1396263};
/// Minimum pitch setting in pitch controller, rad
static constexpr double kPC_MinPit{0.0};
/// Desired (reference) HSS speed for pitch controller, rad/s
static constexpr double kPC_RefSpd{122.9096};
/// Factor to convert radians to degrees
// static constexpr double kR2D{57.295780};
/// Factor to convert radians per second to revolutions per minute
// static constexpr double kRPS2RPM{9.5492966};
/// Rated generator slip percentage in Region 2 1/2, %
static constexpr double kVS_SlPc{10.0};

struct InternalState {
    double generator_speed_filtered;  // Filtered HSS (generator) speed, rad/s.
    double integral_speed_error;      // Current integral of speed error w.r.t. time, rad
    double generator_torque_lastest;  // Commanded electrical generator torque the last time the
                                      // controller was called, N-m
    double time_latest;               // Last time this DLL was called, sec
    double pitch_controller_latest;   // Last time the pitch  controller was called, sec
    double torque_controller_latest;  // Last time the torque controller was called, sec
    std::array<double, 3> pitch_commanded_latest;  // Commanded pitch of each blade the last time the
                                                   // controller was called, rad
    double VS_torque_slope_15;        // Torque/speed slope of region 1 1/2 cut-in torque ramp,
                                      // N-m/(rad/s)
    double VS_torque_slope_25;        // Torque/speed slope of region 2 1/2 induction
                                      // generator, N-m/(rad/s)
    double VS_sync_speed;             // Synchronous speed of region 2 1/2 induction generator, rad/s
    double VS_generator_speed_trans;  // Transitional generator speed (HSS side) between regions
                                      // 2 and 2 1/2, rad/s. 1/2, rad/s
};

inline int SetupFirstCall(const ControllerIO& swap, InternalState& state, char* const avcMSG) {
    // Inform users that we are using this user-defined routine:
    auto aviFAIL = 1;
    strncpy(
        avcMSG,
        "Running with torque and pitch control of the NREL offshore "
        "5MW baseline wind turbine from DISCON.dll as written by J. "
        "Jonkman of NREL/NWTC for use in the IEA Annex XXIII OC3 "
        "studies.",
        swap.message_array_size
    );

    // Determine some torque control parameters not specified directly
    state.VS_sync_speed = kVS_RtGnSp / (1. + 0.01 * kVS_SlPc);
    state.VS_torque_slope_15 = (kVS_Rgn2K * kVS_Rgn2Sp * kVS_Rgn2Sp) / (kVS_Rgn2Sp - kVS_CtInSp);
    state.VS_torque_slope_25 = (kVS_RtPwr / kVS_RtGnSp) / (kVS_RtGnSp - state.VS_sync_speed);
    if (kVS_Rgn2K == 0.) {
        // Region 2 torque is flat, and thus, the denominator in the else condition is zero
        state.VS_generator_speed_trans = state.VS_sync_speed;
    } else {
        // Region 2 torque is quadratic with speed
        state.VS_generator_speed_trans =
            (state.VS_torque_slope_25 -
             sqrt(
                 state.VS_torque_slope_25 *
                 (state.VS_torque_slope_25 - 4.0 * kVS_Rgn2K * state.VS_sync_speed)
             )) /
            (2.0 * kVS_Rgn2K);
    }

    //----------------------------------------------------------------------
    // Check validity of input parameters
    //----------------------------------------------------------------------

    // Initialize aviFAIL to true (will be set to false if all checks pass)
    aviFAIL = -1;

    if (kCornerFreq <= 0.0) {
        strncpy(avcMSG, "CornerFreq must be greater than zero.", swap.message_array_size);
    } else if (kVS_DT <= 0.0) {
        strncpy(avcMSG, "VS_DT must be greater than zero.", swap.message_array_size);
    } else if (kVS_CtInSp < 0.0) {
        strncpy(avcMSG, "VS_CtInSp must not be negative.", swap.message_array_size);
    } else if (kVS_Rgn2Sp <= kVS_CtInSp) {
        strncpy(avcMSG, "VS_Rgn2Sp must be greater than VS_CtInSp.", swap.message_array_size);
    } else if (state.VS_generator_speed_trans < kVS_Rgn2Sp) {
        strncpy(
            avcMSG, "VS_generator_speed_trans must not be less than VS_Rgn2Sp.",
            swap.message_array_size
        );
    } else if (kVS_SlPc <= 0.0) {
        strncpy(avcMSG, "VS_SlPc must be greater than zero.", swap.message_array_size);
    } else if (kVS_MaxRat <= 0.0) {
        strncpy(avcMSG, "VS_MaxRat must be greater than zero.", swap.message_array_size);
    } else if (kVS_RtPwr < 0.0) {
        strncpy(avcMSG, "VS_RtPwr must not be negative.", swap.message_array_size);
    } else if (kVS_Rgn2K < 0.0) {
        strncpy(avcMSG, "VS_Rgn2K must not be negative.", swap.message_array_size);
    } else if (kVS_Rgn2K * kVS_RtGnSp * kVS_RtGnSp > kVS_RtPwr / kVS_RtGnSp) {
        strncpy(
            avcMSG, "VS_Rgn2K*VS_RtGnSp^2 must not be greater than VS_RtPwr/VS_RtGnSp.",
            swap.message_array_size
        );
    } else if (kVS_MaxTq < kVS_RtPwr / kVS_RtGnSp) {
        strncpy(
            avcMSG, "VS_RtPwr/VS_RtGnSp must not be greater than VS_MaxTq.", swap.message_array_size
        );
    } else if (kPC_KI <= 0.0) {
        strncpy(avcMSG, "PC_KI must be greater than zero.", swap.message_array_size);
    } else if (kPC_KK <= 0.0) {
        strncpy(avcMSG, "PC_KK must be greater than zero.", swap.message_array_size);
    } else if (kPC_RefSpd <= 0.0) {
        strncpy(avcMSG, "PC_RefSpd must be greater than zero.", swap.message_array_size);
    } else if (kPC_MaxRat <= 0.0) {
        strncpy(avcMSG, "PC_MaxRat must be greater than zero.", swap.message_array_size);
    } else {
        aviFAIL = 0;
        memset(avcMSG, 0, swap.message_array_size);
    }

    // Initialize the state variables
    // NOTE: generator_torque_lastest is initialized in the torque controller below for
    // simplicity (not here).
    // --------------------------------------------------------------------------------------------
    // This will ensure that generator speed filter will use the initial value of the
    // generator speed on the first pass
    state.generator_speed_filtered = swap.generator_speed_actual;

    // This will ensure that the variable speed controller picks the correct control region
    // and the pitch controller picks the correct gain on the first call
    state.pitch_commanded_latest[0] = swap.pitch_blade1_actual;
    state.pitch_commanded_latest[1] = swap.pitch_blade2_actual;
    state.pitch_commanded_latest[2] = swap.pitch_blade3_actual;

    // This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
    const auto GK = 1.0 / (1.0 + state.pitch_commanded_latest[0] / kPC_KK);

    // This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
    state.integral_speed_error = state.pitch_commanded_latest[1] / (GK * kPC_KI);

    // This will ensure that generator speed filter will use the initial value of the
    // generator speed on the first pass
    state.time_latest = swap.time;

    // This will ensure that the pitch controller is called on the first pass
    state.pitch_controller_latest = swap.time - kPC_DT;

    // This will ensure that the torque controller is called on the first pass
    state.torque_controller_latest = swap.time - kVS_DT;

    return aviFAIL;
}

inline void FilterGeneratorSpeed(const ControllerIO& swap, InternalState& state) {
    // NOTE: This is a very simple recursive, single-pole, low-pass filter with exponential
    // smoothing

    // Update the coefficient in the recursive formula based on the elapsed time since the
    // last call to the controller
    const auto alpha = exp((state.time_latest - swap.time) * kCornerFreq);

    // Apply the filter
    state.generator_speed_filtered =
        (1. - alpha) * swap.generator_speed_actual + alpha * state.generator_speed_filtered;
}

inline void VariableSpeedTorqueControl(ControllerIO& swap, InternalState& state) {
    auto elapsed_time = swap.time - state.torque_controller_latest;

    // Only perform the control calculations if the elapsed time is greater than or equal to
    // the communication interval of the torque controller NOTE: Time is scaled by OnePlusEps
    // to ensure that the controller is called at every time step when kVS_DT = DT, even in
    // the presence of numerical precision errors

    if ((swap.time * kOnePlusEps - state.torque_controller_latest) >= kVS_DT) {
        auto gen_trq = 0.;  // Electrical generator torque, N-m
        // Compute the generator torque, which depends on which region we are in
        if ((state.generator_speed_filtered >= kVS_RtGnSp) ||
            (state.pitch_commanded_latest[0] >= kVS_Rgn3MP)) {
            // We are in region 3 - power is constant
            gen_trq = kVS_RtPwr / state.generator_speed_filtered;
        } else if (state.generator_speed_filtered <= kVS_CtInSp) {
            // We are in region 1 - torque is zero
            gen_trq = 0.0;
        } else if (state.generator_speed_filtered < kVS_Rgn2Sp) {
            // We are in region 1 1/2 - linear ramp in torque from zero to optimal
            gen_trq = state.VS_torque_slope_15 * (state.generator_speed_filtered - kVS_CtInSp);
        } else if (state.generator_speed_filtered < state.VS_generator_speed_trans) {
            // We are in region 2 - optimal torque is proportional to the square of the
            // generator speed
            gen_trq = kVS_Rgn2K * state.generator_speed_filtered * state.generator_speed_filtered;
        } else {
            // We are in region 2 1/2 - simple induction generator transition region
            gen_trq =
                state.VS_torque_slope_25 * (state.generator_speed_filtered - state.VS_sync_speed);
        }

        // Saturate the commanded torque using the maximum torque limit
        gen_trq = std::min(gen_trq, kVS_MaxTq);

        // Initialize the value of generator_torque_lastest on the first pass only
        if (swap.status == 0) {
            state.generator_torque_lastest = gen_trq;
        }

        // Torque rate based on the current and last torque commands, N-m/s
        // Saturate the torque rate using its maximum absolute value
        const auto trq_rate = std::clamp(
            (gen_trq - state.generator_torque_lastest) / elapsed_time, -kVS_MaxRat, kVS_MaxRat
        );

        // Saturate the command using the torque rate limit
        gen_trq = state.generator_torque_lastest + trq_rate * elapsed_time;

        // Reset the values of torque_controller_latest and generator_torque_lastest to the
        // current values
        state.torque_controller_latest = swap.time;
        state.generator_torque_lastest = gen_trq;
    }

    // Set the generator contactor status, avrSWAP(35), to main (high speed)
    //   variable-speed generator, the torque override to yes, and command the
    //   generator torque (See Appendix A of Bladed User's Guide):

    swap.generator_contactor_status =
        1;  // Generator contactor status: 1=main (high speed) variable-speed generator
    swap.torque_override = 0;                                        // Torque override: 0=yes
    swap.generator_torque_command = state.generator_torque_lastest;  // Demanded generator torque
}

inline void PitchControl(ControllerIO& swap, InternalState& state) {
    const auto elapsed_time = swap.time - state.pitch_controller_latest;

    // Only perform the control calculations if the elapsed time is greater than or equal to
    // the communication interval of the pitch controller NOTE: Time is scaled by OnePlusEps
    // to ensure that the contoller is called at every time step when PC_DT = DT, even in the
    // presence of numerical precision errors
    if ((swap.time * kOnePlusEps - state.pitch_controller_latest) >= kPC_DT) {
        // Current value of the gain correction factor, used in the gain
        // scheduling law of the pitch controller, (-).
        // Based on the previously commanded pitch angle for blade 1:
        const auto GK = 1.0 / (1.0 + state.pitch_commanded_latest[0] / kPC_KK);

        // Compute the current speed error and its integral w.r.t. time; saturate the
        // integral term using the pitch angle limits
        const auto speed_error = state.generator_speed_filtered - kPC_RefSpd;
        state.integral_speed_error += speed_error * elapsed_time;
        state.integral_speed_error = std::clamp(
            state.integral_speed_error, kPC_MinPit / (kOnePlusEps * kPC_KI),
            kPC_MaxPit / (kOnePlusEps * kPC_KI)
        );

        // Compute the pitch commands associated with the proportional and integral gains
        const auto pitch_com_proportional = GK * kPC_KP * speed_error;
        const auto pitch_com_integral = GK * kPC_KI * state.integral_speed_error;

        // Superimpose the individual commands to get the total pitch command; saturate the
        // overall command using the pitch angle limits
        const auto pitch_com_total =
            std::clamp(pitch_com_proportional + pitch_com_integral, kPC_MinPit, kPC_MaxPit);
        // Saturate the overall commanded pitch using the pitch rate limit:
        // NOTE: Since the current pitch angle may be different for each blade
        //       (depending on the type of actuator implemented in the structural
        //       dynamics model), this pitch rate limit calculation and the
        //       resulting overall pitch angle command may be different for each
        //       blade.

        // Current values of the blade pitch angles, rad
        const auto blade_pitch =
            std::array{swap.pitch_blade1_actual, swap.pitch_blade2_actual, swap.pitch_blade3_actual};

        // Pitch rates of each blade based on the current pitch angles and current pitch
        // command, rad/s
        auto pitch_rate = std::array<double, 3>{};

        // Loop through all blades
        std::ranges::transform(blade_pitch, std::begin(pitch_rate), [&](auto pitch) {
            auto rate = (pitch_com_total - pitch) / elapsed_time;
            return std::clamp(rate, -kPC_MaxRat, kPC_MaxRat);
        });

        std::ranges::transform(
            blade_pitch, pitch_rate, std::begin(state.pitch_commanded_latest),
            [&](auto pitch, auto rate) {
                auto commanded = pitch + (rate * elapsed_time);
                return std::clamp(commanded, kPC_MinPit, kPC_MaxPit);
            }
        );

        // Reset the value of pitch_controller_latest to the current value
        state.pitch_controller_latest = swap.time;
    }

    // Set the pitch override to yes and command the pitch demanded from the last
    // call to the controller (See Appendix A of Bladed User's Guide):
    swap.pitch_override = 0.;  // Pitch override: 0=yes

    swap.pitch_blade1_command = state.pitch_commanded_latest[0];  // Use the command angles of all
                                                                  // blades if using individual pitch
    swap.pitch_blade2_command = state.pitch_commanded_latest[1];  // "
    swap.pitch_blade3_command = state.pitch_commanded_latest[2];  // "

    swap.pitch_collective_command =
        state.pitch_commanded_latest[0];  // Use the command angle of blade 1 if using collective
                                          // pitch
}

inline int ComputeControl(ControllerIO& swap, InternalState& state, char* const avcMSG) {
    // Abort if the user has not requested a pitch angle actuator (See Appendix A of Bladed
    // User's Guide)
    auto aviFAIL = 0;
    if (swap.pitch_actuator_type_req != 0) {
        aviFAIL = -1;
        strncpy(avcMSG, "Pitch angle actuator not requested.", swap.message_array_size);
    }

    // Set unused outputs to zero (See Appendix A of Bladed User's Guide):
    swap.shaft_brake_status = 0;             // Shaft brake status: 0=off
    swap.yaw_actuator_torque_command = 0.;   // Demanded yaw actuator torque
    swap.pitch_rate_command = 0.;            // Demanded pitch rate (Collective pitch)
    swap.nacelle_yaw_rate_command = 0.;      // Demanded nacelle yaw rate
    swap.n_log_variables = 0U;               // Number of variables returned for logging
    swap.generator_startup_resistance = 0.;  // Generator start-up resistance
    swap.loads_request = 0;                  // Request for loads: 0=none
    swap.variable_slip_status = 0;           // Variable slip current status
    swap.variable_slip_demand = 0;           // Variable slip current demand

    FilterGeneratorSpeed(swap, state);

    VariableSpeedTorqueControl(swap, state);

    PitchControl(swap, state);

    state.time_latest = swap.time;

    return aviFAIL;
}

inline void PackInternalStateToFile(const InternalState& state, const char* const accINFILE) {
    auto fp = std::ofstream(accINFILE, std::ios::binary);
    if (fp) {
        fp << state.generator_speed_filtered << state.integral_speed_error
           << state.generator_torque_lastest << state.time_latest << state.pitch_controller_latest
           << state.torque_controller_latest << state.pitch_commanded_latest[0]
           << state.pitch_commanded_latest[1] << state.pitch_commanded_latest[2]
           << state.VS_torque_slope_15 << state.VS_torque_slope_25 << state.VS_sync_speed
           << state.VS_generator_speed_trans;
    }
}

inline void UnpackInternalStateFromFile(const char* const accINFILE, InternalState& state) {
    auto fp = std::ifstream(accINFILE, std::ios::binary);
    if (fp) {
        fp >> state.generator_speed_filtered;
        fp >> state.integral_speed_error;
        fp >> state.generator_torque_lastest;
        fp >> state.time_latest;
        fp >> state.pitch_controller_latest;
        fp >> state.torque_controller_latest;
        fp >> state.pitch_commanded_latest[0];
        fp >> state.pitch_commanded_latest[1];
        fp >> state.pitch_commanded_latest[2];
        fp >> state.VS_torque_slope_15;
        fp >> state.VS_torque_slope_25;
        fp >> state.VS_sync_speed;
        fp >> state.VS_generator_speed_trans;
    }
}
// TODO This is a quick and dirty conversion of the DISCON function from the original C code to
// C++. It needs to be refactored to be more idiomatic C++.
void DISCON(
    float avrSWAP[], int* aviFAIL, const char* const accINFILE, char* const, char* const avcMSG
) {
    // Internal state
    static InternalState state;

    // Map swap from calling program to struct
    ControllerIO swap;
    auto swap_array = std::array<float, kSwapArraySize>{};
    std::copy(avrSWAP, std::next(avrSWAP, 81), swap_array.begin());
    swap.CopyFromSwapArray(swap_array);

    // Initialize aviFAIL to 0
    *aviFAIL = 0;

    if (swap.status == 0) {
        *aviFAIL = SetupFirstCall(swap, state, avcMSG);
    }

    if (swap.status >= 0 && *aviFAIL >= 0) {
        *aviFAIL = ComputeControl(swap, state, avcMSG);
    } else if (swap.status == -8) {
        PackInternalStateToFile(state, accINFILE);
    } else if (swap.status == -9) {
        UnpackInternalStateFromFile(accINFILE, state);
    }

    swap.CopyToSwapArray(swap_array);
    std::ranges::copy(swap_array, avrSWAP);
}

}  // extern "C"

}  // namespace kynema::util
