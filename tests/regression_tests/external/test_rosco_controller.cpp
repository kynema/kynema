#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "interfaces/components/controller.hpp"

#include "Kynema_config.h"

namespace kynema::tests {

TEST(ROSCO_Controller, initialize) {
    const auto shared_lib_path = std::string{static_cast<const char*>(Kynema_ROSCO_LIBRARY)};
    const auto controller_function_name = std::string{"DISCON"};

    auto controller = interfaces::components::Controller(interfaces::components::ControllerInput{
        .shared_lib_path = shared_lib_path,
        .function_name = controller_function_name,
        .input_file_path = "./IEA-15-240-RWT/DISCON.IN",
        .output_file_path = ""
    });

    controller.io.status = 0;
    controller.io.time = 0.;
    controller.io.dt = 0.01;
    controller.io.rotor_speed_actual = 5.;

    controller.CallController();
}

}  // namespace kynema::tests
