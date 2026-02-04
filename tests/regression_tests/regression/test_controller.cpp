#include <array>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "interfaces/components/controller.hpp"
#include "interfaces/components/controller_input.hpp"
#include "interfaces/components/controller_io.hpp"
#include "vendor/dylib/dylib.hpp"

namespace kynema::tests {

TEST(ControllerTest, DisconController) {
    // Test data generated using the following regression test from the
    // OpenFAST/r-test repository:
    // https://github.com/OpenFAST/r-test/tree/main/glue-codes/openfast/5MW_Land_DLL_WTurb
    // at time = 0.0s

    // Use dylib to load the dynamic library and get access to the controller functions
    const util::dylib lib("./DISCON.dll", util::dylib::no_filename_decorations);
    auto DISCON = lib.get_function<void(float*, int&, char*, char*, char*)>("DISCON");

    interfaces::components::ControllerIO swap;
    swap.status = 0.;
    swap.time = 0.;
    swap.pitch_blade1_actual = 0.;
    swap.pitch_actuator_type_req = 0.;
    swap.generator_speed_actual = 122.909576;
    swap.horizontal_wind_speed = 11.9900799;
    swap.pitch_blade2_actual = 0.;
    swap.pitch_blade3_actual = 0.;
    swap.generator_contactor_status = 1.;
    swap.shaft_brake_status = 0.;
    swap.yaw_actuator_torque_command = 0.;
    swap.pitch_blade1_command = 0.;
    swap.pitch_blade2_command = 0.;
    swap.pitch_blade3_command = 0.;
    swap.pitch_collective_command = 0.;
    swap.pitch_rate_command = 0.;
    swap.generator_torque_command = 0.;
    swap.nacelle_yaw_rate_command = 0.;
    swap.message_array_size = 3.;
    swap.infile_array_size = 82.;
    swap.outname_array_size = 96.;
    swap.pitch_override = 0.;
    swap.torque_override = 0.;
    swap.n_blades = 3.;
    swap.n_log_variables = 0.;
    swap.generator_startup_resistance = 0.;
    swap.loads_request = 0.;
    swap.variable_slip_status = 0.;
    swap.variable_slip_demand = 0.;

    auto avrSWAP = std::array<float, interfaces::components::kSwapArraySize>{};
    swap.CopyToSwapArray(avrSWAP);

    // Expect demanded generator torque to be 0. before calling the controller
    EXPECT_FLOAT_EQ(avrSWAP[47], 0.);

    // Call DISCON and expect the following outputs
    int aviFAIL = 0;
    char in_file[] = "in_file";
    char out_name[] = "out_name";
    char msg[] = "msg";
    DISCON(
        avrSWAP.data(), aviFAIL, static_cast<char*>(in_file), static_cast<char*>(out_name),
        static_cast<char*>(msg)
    );

    EXPECT_FLOAT_EQ(avrSWAP[34], 1.);           // GeneratorContactorStatus
    EXPECT_FLOAT_EQ(avrSWAP[35], 0.);           // ShaftBrakeStatus
    EXPECT_FLOAT_EQ(avrSWAP[40], 0.);           // DemandedYawActuatorTorque
    EXPECT_FLOAT_EQ(avrSWAP[44], 0.);           // PitchComCol
    EXPECT_FLOAT_EQ(avrSWAP[46], 43093.5508F);  // DemandedGeneratorTorque
    EXPECT_FLOAT_EQ(avrSWAP[47], 0.);           // DemandedNacelleYawRate
}

TEST(ControllerTest, Controller) {
    // Get a handle to the controller function via the Controller class and use it to
    // calculate the controller outputs
    const auto shared_lib_path = std::string{"./DISCON.dll"};
    const auto controller_function_name = std::string{"DISCON"};

    auto controller = interfaces::components::Controller(interfaces::components::ControllerInput{
        .shared_lib_path = shared_lib_path,
        .function_name = controller_function_name,
        .input_file_path = "",
        .output_file_path = ""
    });

    controller.io.status = 0.;
    controller.io.time = 0.;
    controller.io.pitch_blade1_actual = 0.;
    controller.io.pitch_actuator_type_req = 0.;
    controller.io.generator_speed_actual = 122.909576;
    controller.io.horizontal_wind_speed = 11.9900799;
    controller.io.pitch_blade2_actual = 0.;
    controller.io.pitch_blade3_actual = 0.;
    controller.io.generator_contactor_status = 1.;
    controller.io.shaft_brake_status = 0.;
    controller.io.yaw_actuator_torque_command = 0.;
    controller.io.pitch_blade1_command = 0.;
    controller.io.pitch_blade2_command = 0.;
    controller.io.pitch_blade3_command = 0.;
    controller.io.pitch_collective_command = 0.;
    controller.io.pitch_rate_command = 0.;
    controller.io.generator_torque_command = 0.;
    controller.io.nacelle_yaw_rate_command = 0.;
    controller.io.message_array_size = 3.;
    controller.io.infile_array_size = 82.;
    controller.io.outname_array_size = 96.;
    controller.io.pitch_override = 0.;
    controller.io.torque_override = 0.;
    controller.io.n_blades = 3.;
    controller.io.n_log_variables = 0.;
    controller.io.generator_startup_resistance = 0.;
    controller.io.loads_request = 0.;
    controller.io.variable_slip_status = 0.;
    controller.io.variable_slip_demand = 0.;

    EXPECT_DOUBLE_EQ(controller.io.generator_torque_command, 0.);

    controller.CallController();

    EXPECT_DOUBLE_EQ(controller.io.generator_contactor_status, 1.);   // GeneratorContactorStatus
    EXPECT_DOUBLE_EQ(controller.io.shaft_brake_status, 0.);           // ShaftBrakeStatus
    EXPECT_DOUBLE_EQ(controller.io.yaw_actuator_torque_command, 0.);  // DemandedYawActuatorTorque
    EXPECT_DOUBLE_EQ(controller.io.pitch_collective_command, 0.);     // PitchComCol
    EXPECT_DOUBLE_EQ(
        controller.io.generator_torque_command, 43093.55078125
    );                                                             // DemandedGeneratorTorque
    EXPECT_DOUBLE_EQ(controller.io.nacelle_yaw_rate_command, 0.);  // DemandedNacelleYawRate
}

TEST(ControllerTest, ControllerExceptionInvalidSharedLibraryPath) {
    // Test case: invalid shared library path
    const auto shared_lib_path = std::string{"./INVALID.dll"};
    const auto controller_function_name = std::string{"DISCON"};

    EXPECT_THROW(
        auto controller = interfaces::components::Controller(interfaces::components::ControllerInput{
            .shared_lib_path = shared_lib_path,
            .function_name = controller_function_name,
            .input_file_path = "",
            .output_file_path = ""
        }),
        std::runtime_error
    );
}

TEST(ControllerTest, ControllerExceptionInvalidControllerFunctionName) {
    // Test case: invalid controller function name
    const auto shared_lib_path = std::string{"./DISCON.dll"};
    const auto controller_function_name = std::string{"INVALID"};

    EXPECT_THROW(
        auto controller = interfaces::components::Controller(interfaces::components::ControllerInput{
            .shared_lib_path = shared_lib_path,
            .function_name = controller_function_name,
            .input_file_path = "",
            .output_file_path = ""
        }),
        std::runtime_error
    );
}

}  // namespace kynema::tests
