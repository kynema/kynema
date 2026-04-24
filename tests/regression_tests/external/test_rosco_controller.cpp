#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "utilities/controllers/turbine_controller.hpp"

#include "Kynema_config.h"

namespace kynema_fmb::tests {

TEST(ROSCO_Controller, initialize) {
    const auto shared_lib_path = std::string{static_cast<const char*>(Kynema_ROSCO_LIBRARY)};
    const auto controller_function_name = std::string{"DISCON"};

    auto controller = util::TurbineController(
        shared_lib_path, controller_function_name, "./IEA-15-240-RWT/DISCON.IN", ""
    );

    controller.io.status = 0;
    controller.io.time = 0.;
    controller.io.dt = 0.01;
    controller.io.rotor_speed_actual = 5.;

    controller.CallController();
}

}  // namespace kynema_fmb::tests
