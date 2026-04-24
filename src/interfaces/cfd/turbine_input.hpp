#pragma once

#include "interfaces/cfd/floating_platform_input.hpp"

namespace kynema_fmb::interfaces::cfd {

/**
 * @brief A collection of the input objects defining the CFD problem's configuration
 */
struct TurbineInput {
    // Floating platform
    FloatingPlatformInput floating_platform;
};

}  // namespace kynema_fmb::interfaces::cfd
