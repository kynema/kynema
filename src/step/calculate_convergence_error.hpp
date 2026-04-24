#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "constraints/constraints.hpp"
#include "solver/calculate_error_sum_squares.hpp"
#include "solver/solver.hpp"
#include "state/state.hpp"
#include "step/step_parameters.hpp"

namespace kynema_fmb::step {

/// @brief Calculation based on Table 1 of DOI: 10.1115/1.4033441
template <typename DeviceType>
inline double CalculateConvergenceError(
    const StepParameters& parameters, const Solver<DeviceType>& solver,
    const State<DeviceType>& state, const Constraints<DeviceType>& constraints
) {
    auto region = Kokkos::Profiling::ScopedRegion("Calculate Convergence Error");

    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;

    auto sum_error_squared_system = 0.;
    Kokkos::parallel_reduce(
        RangePolicy(0, solver.num_system_nodes),
        solver::CalculateSystemErrorSumSquares<DeviceType>{
            parameters.absolute_convergence_tol,
            parameters.relative_convergence_tol,
            parameters.h,
            state.active_dofs,
            state.node_freedom_map_table,
            state.q_delta,
            solver.x,
        },
        sum_error_squared_system
    );
    auto sum_error_squared_constraints = 0.;
    Kokkos::parallel_reduce(
        RangePolicy(0, constraints.num_constraints),
        solver::CalculateConstraintsErrorSumSquares<DeviceType>{
            parameters.absolute_convergence_tol,
            parameters.relative_convergence_tol,
            solver.num_system_dofs,
            constraints.row_range,
            constraints.lambda,
            solver.x,
        },
        sum_error_squared_constraints
    );

    const auto sum_error_squared = sum_error_squared_system + sum_error_squared_constraints;
    return std::sqrt(sum_error_squared / static_cast<double>(solver.num_dofs));
}

}  // namespace kynema_fmb::step
