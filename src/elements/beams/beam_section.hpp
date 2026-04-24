#pragma once

#include <array>

namespace kynema_fmb {

/**
 * @brief A structure containing the position, mass matrix, and stiffness matrix to be used
 * in defining a beam cross section at a given location.
 */
struct BeamSection {
    double position;                              // Position of section in element on range [0, 1]
    std::array<std::array<double, 6>, 6> M_star;  // Mass matrix in material frame
    std::array<std::array<double, 6>, 6> C_star;  // Stiffness matrix in material frame

    BeamSection(
        double s, const std::array<std::array<double, 6>, 6>& mass,
        const std::array<std::array<double, 6>, 6>& stiffness
    )
        : position(s), M_star(mass), C_star(stiffness) {}
};

}  // namespace kynema_fmb
