#pragma once

#include "solution_input.hpp"

namespace kynema::interfaces::components {

/**
 * @brief A factory class for building a solution input, which is used by interfaces to control
 * the time stepping proceedure.
 *
 * @details Each of the methods returns a reference to this SolutionBuilder object, which allows
 * chaining multiple calls together in a single statement.
 */
struct SolutionBuilder {
    /**
     * @brief solves this problem with the static solver
     *
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& EnableStaticSolve();

    /**
     * @brief Solves this problem with the dynamic solver
     *
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& EnableDynamicSolve();

    /**
     * @brief Sets the gravity vector for the problem
     *
     * @param gravity The gravity vector
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& SetGravity(const std::array<double, 3>& gravity);

    /**
     * @brief Sets the relative error used to determine if nonlinear iterations have converged
     *
     * @param rtol the relative error tolerance
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& SetRelativeErrorTolerance(double rtol);

    /**
     * @brief Sets the absolute error used to determine if nonlinear iterations have converged
     *
     * @param atol The absolute error tolerance
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& SetAbsoluteErrorTolerance(double atol);

    /**
     * @brief Sets the maximum number of nonlinear iterations to take at each time step
     *
     * @param max_iter The maximum number of nonlinear iterations
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& SetMaximumNonlinearIterations(size_t max_iter);

    /**
     * @brief Sets the timestep size
     *
     * @param time_step The time step size
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& SetTimeStep(double time_step);

    /**
     * @brief Sets the numerical damping factor for the generalized alpha solver
     *
     * @details Ranged from 0 to 1
     *
     * @param rho_inf The damping factor
     * @return A reference to this solution builder object to allow chaining
     */
    SolutionBuilder& SetDampingFactor(double rho_inf);

    /**
     * @brief Creates a SolutionInput object based on the previously set parameters
     *
     * @return The completed SolutionInput object
     */
    [[nodiscard]] const SolutionInput& Input() const;

private:
    SolutionInput input;
};

}  // namespace kynema::interfaces::components
