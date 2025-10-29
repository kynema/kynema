#pragma once

#include "interfaces/components/aerodynamics_builder.hpp"
#include "interfaces/components/controller_builder.hpp"
#include "interfaces/components/outputs_builder.hpp"
#include "interfaces/components/solution_builder.hpp"
#include "interfaces/components/turbine_builder.hpp"
#include "turbine_interface.hpp"

namespace kynema::interfaces {

/**
 * @brief Builder class to construct a TurbineInterface by composing Turbine and Solution components
 *
 * @details This class combines these builders through a facade pattern, providing a unified
 * interface for constructing a complete blade model while keeping the configuration of
 * individual components separate and maintainable.
 *
 * - SolutionBuilder: Configures solver type, tolerances, time steps, and output settings
 * - TurbineBuilder: Configures blade geometry, reference axes, section properties, and structural
 *                 matrices
 */
class TurbineInterfaceBuilder {
public:
    /// @brief Gets the builder for the solution component
    /// @return A reference to the SolutionBuilder for the solution component
    [[nodiscard]] components::SolutionBuilder& Solution() { return this->solution_builder; }

    /// @brief Get the builder for the turbine component
    /// @return A reference to the TurbineBuilder for the turbine component
    [[nodiscard]] components::TurbineBuilder& Turbine() { return this->turbine_builder; }

    /// @brief Get the builder for the aerodynamics component
    /// @return a reference to the AerodynamicsBuilder for the aerodynamics component
    [[nodiscard]] components::AerodynamicsBuilder& Aerodynamics() { return this->aero_builder; }

    /// @brief Get the builder for the controller component
    /// @return a reference to the ControllerBuilder for the controller component
    [[nodiscard]] components::ControllerBuilder& Controller() { return this->controller_builder; }

    /// @brief Gets the builder for the outputs component
    /// @return A reference to the OutputsBuilder attribute
    [[nodiscard]] components::OutputsBuilder& Outputs() { return this->outputs_builder; }

    /**
     * @brief Builds the TurbineInterface by composing the blade, tower, nacelle, hub, and solution
     * components
     * @return A TurbineInterface object
     */
    [[nodiscard]] TurbineInterface Build() {
        return TurbineInterface(
            this->solution_builder.Input(), this->turbine_builder.Input(),
            this->aero_builder.Input(), this->controller_builder.Input(),
            this->outputs_builder.Config()
        );
    }

private:
    components::SolutionBuilder solution_builder;      ///< Builder for the Solution component
    components::TurbineBuilder turbine_builder;        ///< Builder for the Turbine component
    components::AerodynamicsBuilder aero_builder;      ///< Builder for the Aerodynamics component
    components::ControllerBuilder controller_builder;  ///< Builder for the Controller component
    components::OutputsBuilder outputs_builder;        ///< Builder for the Outputs component
};

}  // namespace kynema::interfaces
