#include "discon_rotor_test_controller.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iterator>
#include <numbers>

#include "interfaces/components/controller_io.hpp"

namespace kynema::util {

extern "C" {

void PITCH_CONTROLLER(
    float avrSWAP[], int* aviFAIL, const char* const, char* const, char* const avcMSG
) {
    static auto first_call = true;
    // Map swap from calling program to struct
    interfaces::components::ControllerIO io;
    auto swap_array = std::array<float, interfaces::components::kSwapArraySize>{};
    std::copy(avrSWAP, std::next(avrSWAP, 81), swap_array.begin());
    io.CopyFromSwapArray(swap_array);

    // Update pitch angle in radians
    io.pitch_blade1_command = 2. * std::numbers::pi * io.time / 3.;  // Full rotation every 3 seconds
    io.pitch_blade2_command = 2. * std::numbers::pi * io.time / 6.;  // Full rotation every 6 seconds
    io.pitch_blade3_command = 2. * std::numbers::pi * io.time / 9.;  // Full rotation every 9 seconds

    // Set failure flag to zero
    *aviFAIL = 0;

    // Set message to success
    strncpy(avcMSG, "success\0", io.message_array_size);
    *std::next(
        avcMSG,
        static_cast<typename std::iterator_traits<char*>::difference_type>(io.message_array_size) - 1
    ) = 0;

    first_call = false;

    io.CopyToSwapArray(swap_array);
    std::ranges::copy(swap_array, avrSWAP);
}

}  // extern "C"

}  // namespace kynema::util
