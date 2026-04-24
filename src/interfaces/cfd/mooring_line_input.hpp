#pragma once

#include <array>

namespace kynema_fmb::interfaces::cfd {

/**
 * @brief A descritpion of the configuration of a mooring line for use in initialization
 */
struct MooringLineInput {
    /// Mooring line stiffness
    double stiffness = 0.;

    /// Undeformed length of mooring line
    double undeformed_length = 0.;

    /// Fairlead point coordinates (XYZ)
    std::array<double, 3> fairlead_position{0., 0., 0.};

    /// Fairlead point velocity (XYZ)
    std::array<double, 3> fairlead_velocity{0., 0., 0.};

    /// Fairlead point acceleration (XYZ)
    std::array<double, 3> fairlead_acceleration{0., 0., 0.};

    /// Anchor point coordinates (XYZ)
    std::array<double, 3> anchor_position{0., 0., 0.};

    /// Anchor point velocity (XYZ)
    std::array<double, 3> anchor_velocity{0., 0., 0.};

    /// Anchor point acceleration (XYZ)
    std::array<double, 3> anchor_acceleration{0., 0., 0.};
};

}  // namespace kynema_fmb::interfaces::cfd
