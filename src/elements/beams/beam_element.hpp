#pragma once

#include <array>
#include <span>
#include <vector>

#include "beam_section.hpp"

namespace kynema_fmb {

/**
 * @brief Beam element constitutes flexible beams material behavior in kynema.
 *
 * @details A beam element is defined by a set of nodes and sections. Each section is defined by a
 * 6x6 mass matrix and a 6x6 stiffness matrix. The element quadrature is used to integrate the
 * mass and stiffness matrices along the length of the beam.
 */
struct BeamElement {
    size_t ID;                                      // Element identifier
    std::vector<size_t> node_ids;                   // Element node identifiers
    std::vector<BeamSection> sections;              // Element mass/stiffness in material frame
    std::vector<std::array<double, 2>> quadrature;  // Element quadrature points and weights
    std::array<double, 6> mu;  // Element stiffness-proportional damping coefficients

    BeamElement(
        size_t id, std::span<const size_t> n, std::span<const BeamSection> s,
        std::span<const std::array<double, 2>> q,
        const std::array<double, 6>& mu_ = {0., 0., 0., 0., 0., 0.}
    )
        : ID(id), mu(mu_) {
        node_ids.assign(std::begin(n), std::end(n));
        sections.assign(std::begin(s), std::end(s));
        quadrature.assign(std::begin(q), std::end(q));
    }
};

}  // namespace kynema_fmb
