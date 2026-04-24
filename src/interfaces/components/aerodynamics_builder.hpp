#pragma once

#include <ranges>
#include <span>

#include "aerodynamics_input.hpp"

namespace kynema_fmb::interfaces::components {

class AerodynamicsBuilder {
public:
    AerodynamicsBuilder() = default;

    AerodynamicsInput& Input() { return input; }

    AerodynamicsBuilder& EnableAero() {
        input.is_enabled = true;
        return *this;
    }

    AerodynamicsBuilder& SetNumberOfAirfoils(size_t number_of_blades) {
        input.aero_inputs.resize(number_of_blades);
        return *this;
    }

    AerodynamicsBuilder& SetAirfoilToBladeMap(std::span<const size_t> map) {
        input.airfoil_map.resize(map.size());
        std::ranges::copy(map, std::begin(input.airfoil_map));
        return *this;
    }

    AerodynamicsBuilder& SetAirfoilSections(
        size_t airfoil_number, std::span<const AerodynamicSection> sections
    ) {
        input.aero_inputs[airfoil_number].resize(sections.size());
        std::ranges::copy(sections, std::begin(input.aero_inputs[airfoil_number]));
        return *this;
    }

private:
    AerodynamicsInput input;
};

}  // namespace kynema_fmb::interfaces::components
