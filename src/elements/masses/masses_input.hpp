#pragma once

#include <array>
#include <span>
#include <vector>

#include "mass_element.hpp"

namespace kynema_fmb {

/**
 * @brief Represents the input data for creating mass/rigid body elements
 */
struct MassesInput {
    std::vector<MassElement> elements;  //< All elements present in the masses/rigid_body
    std::array<double, 3> gravity;      //< Gravity vector

    MassesInput(std::span<const MassElement> elems, std::span<const double, 3> g)
        : gravity({g[0], g[1], g[2]}) {
        elements.assign(std::begin(elems), std::end(elems));
    }

    /// Returns the total number of elements present in masses/rigid_body portion of the mesh
    [[nodiscard]] size_t NumElements() const { return elements.size(); }
};

}  // namespace kynema_fmb
