#pragma once

#include "interfaces/components/beam_builder.hpp"
#include "interfaces/components/outputs_builder.hpp"
#include "interfaces/components/solution_builder.hpp"

namespace kynema::interfaces {

class BladeInterface;

/**
 * @brief Builder class to construct a BladeInterface by composing Blade and Solution components
 *
 * @details This class combines these builders through a facade pattern, providing a unified
 * interface for constructing a complete blade model while keeping the configuration of
 * individual components separate and maintainable.
 *
 * - SolutionBuilder: Configures solver type, tolerances, time steps, and output settings
 * - BladeBuilder: Configures blade geometry, reference axes, section properties, and structural
 *                 matrices
 */
class BladeInterfaceBuilder {
public:
    /**
     * @brief Returns the owned SolutionBuilder for defining the solution control variables
     *
     * @return A reference to the SolutionBuilder object
     */
    [[nodiscard]] components::SolutionBuilder& Solution();

    /**
     * @brief Returns the owned BeamBuilder for use in defining the blade object
     *
     * @return A reference to the BeamBuilder object
     */
    [[nodiscard]] components::BeamBuilder& Blade();

    /**
     * @brief Returns the owned OutputsBuilder for use in defining the outputs configuration
     *
     * @return A reference to the OutputsBuilder object
     */
    [[nodiscard]] components::OutputsBuilder& Outputs();

    /**
     * @brief Builds the BladeInterface by composing the Blade and Solution components
     * @return A BladeInterface object
     */
    [[nodiscard]] BladeInterface Build() const;

private:
    components::SolutionBuilder solution_builder;  ///< Builder for the Solution component
    components::BeamBuilder beam_builder;          ///< Builder for the Blade component
    components::OutputsBuilder outputs_builder;    ///< Builder for the Outputs component
};

}  // namespace kynema::interfaces
