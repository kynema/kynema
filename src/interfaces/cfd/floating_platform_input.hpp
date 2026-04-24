#pragma once

#include <array>
#include <vector>

#include "interfaces/cfd/mooring_line_input.hpp"

namespace kynema_fmb::interfaces::cfd {

/**
 * @brief The input configuration options describing a FloatingPlatform object
 */
struct FloatingPlatformInput {
    /// @brief Flag to enable use of floating platform in model
    bool enable{false};

    /// @brief Platform point coordinates and orientation (XYZ,quaternion)
    std::array<double, 7> position{0., 0., 0., 1., 0., 0., 0.};

    /// @brief Platform point translational and rotational velocity
    std::array<double, 6> velocity{0., 0., 0., 0., 0., 0.};

    /// @brief Platform point translational and rotational acceleration
    std::array<double, 6> acceleration{0., 0., 0., 0., 0., 0.};

    /// @brief Platform point mass matrix
    std::array<std::array<double, 6>, 6> mass_matrix{};

    /// @brief Mooring line array
    std::vector<MooringLineInput> mooring_lines;
};

}  // namespace kynema_fmb::interfaces::cfd
