#pragma once

#include <array>
#include <string>

#include "step/step_parameters.hpp"

namespace kynema::interfaces::components {

/**
 * @brief A configuration object used to create the low level StepParameters object
 */
struct SolutionInput {
    /// @brief Array of gravity components (XYZ)
    std::array<double, 3> gravity{0., 0., 0.};

    /// @brief Flag to toggle between static and dynamic solve
    bool dynamic_solve{true};

    /// @brief Solver time step
    double time_step{0.01};

    /// @brief Solver numerical damping factor (0 = maximum damping)
    double rho_inf{0.0};

    /// @brief Maximum number of convergence iterations
    size_t max_iter{5U};

    /// @brief Absolute error tolerance
    double absolute_error_tolerance{1e-5};

    /// @brief Relative error tolerance
    double relative_error_tolerance{1e-3};

    /// @brief  Construct step parameters from inputs
    /// @return Step parameters struct
    [[nodiscard]] StepParameters Parameters() const {
        return {
            this->dynamic_solve,
            this->max_iter,
            this->time_step,
            this->rho_inf,
            this->absolute_error_tolerance,
            this->relative_error_tolerance
        };
    }
};

}  // namespace kynema::interfaces::components
