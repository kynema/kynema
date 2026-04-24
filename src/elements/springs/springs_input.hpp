#pragma once

#include <span>
#include <vector>

#include "spring_element.hpp"

namespace kynema_fmb {

/**
 * @brief Represents the input data for creating spring elements
 */
struct SpringsInput {
    std::vector<SpringElement> elements;  //< All spring elements in the system

    explicit SpringsInput(std::span<const SpringElement> elems) {
        elements.assign(std::begin(elems), std::end(elems));
    }

    /// Returns the total number of spring elements in the system
    [[nodiscard]] size_t NumElements() const { return elements.size(); }
};

}  // namespace kynema_fmb
