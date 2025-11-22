#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "assemble_constraints_matrix.hpp"
#include "assemble_constraints_residual.hpp"
#include "assemble_system_matrix.hpp"
#include "assemble_system_residual.hpp"
#include "calculate_convergence_error.hpp"
#include "constraints/calculate_constraint_output.hpp"
#include "constraints/constraints.hpp"
#include "elements/elements.hpp"
#include "predict_next_state.hpp"
#include "reset_constraints.hpp"
#include "reset_solver.hpp"
#include "solve_system.hpp"
#include "solver/solver.hpp"
#include "state/clone_state.hpp"
#include "state/copy_state_data.hpp"
#include "state/state.hpp"
#include "state/update_algorithmic_acceleration.hpp"
#include "state/update_global_position.hpp"
#include "step_parameters.hpp"
#include "update_constraint_prediction.hpp"
#include "update_constraint_variables.hpp"
#include "update_state_prediction.hpp"
#include "update_system_variables.hpp"
#include "update_tangent_operator.hpp"

namespace kynema {

/**
 * @brief Attempts to complete a single time step in the dynamic/static FEA simulation
 *
 * @param parameters Simulation step parameters including time step size and convergence criteria
 * @param solver     Solver object containing system matrices and solution methods
 * @param elements   Collection of elements (beams, masses etc.) in the FE mesh
 * @param state      Current state of the system (positions, velocities, accelerations etc.)
 * @param constraints System constraints and their associated data
 *
 * @return true if the step converged within the maximum allowed iterations, otherwise false
 */
template <typename DeviceType>
inline bool SolveStep(
    StepParameters& parameters, Solver<DeviceType>& solver, Elements<DeviceType>& elements,
    State<DeviceType>& state, Constraints<DeviceType>& constraints
) {
    auto region = Kokkos::Profiling::ScopedRegion("SolveStep");

    step::PredictNextState(parameters, state);
    step::ResetConstraints(constraints);

    solver.convergence_err.clear();
    double err{1000.};
    for (auto iter = 0U; err > 1.; ++iter) {
        if (iter >= parameters.max_iter) {
            return false;
        }
        step::ResetSolver(solver);

        step::UpdateTangentOperator(parameters, state);

        step::UpdateSystemVariables(parameters, elements, state);

        step::AssembleSystemResidual(solver, elements, state);

        step::AssembleSystemMatrix(parameters, solver, elements);

        step::UpdateConstraintVariables(state, constraints);

        step::AssembleConstraintsMatrix(solver, constraints);

        step::AssembleConstraintsResidual(solver, constraints);

        step::SolveSystem(parameters, solver);

        err = step::CalculateConvergenceError(parameters, solver, state, constraints);

        solver.convergence_err.push_back(err);

        step::UpdateStatePrediction(parameters, solver, state);

        step::UpdateConstraintPrediction(solver, constraints);
    }

    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;

    auto system_range = RangePolicy(0, solver.num_system_nodes);
    Kokkos::parallel_for(
        "UpdateAlgorithmicAcceleration", system_range,
        state::UpdateAlgorithmicAcceleration<DeviceType>{
            state.a,
            state.vd,
            parameters.alpha_f,
            parameters.alpha_m,
        }
    );

    Kokkos::parallel_for(
        "UpdateGlobalPosition", system_range,
        state::UpdateGlobalPosition<DeviceType>{
            state.q,
            state.x0,
            state.x,
        }
    );

    auto constraints_range = RangePolicy(0, constraints.num_constraints);
    Kokkos::parallel_for(
        "CalculateConstraintOutput", constraints_range,
        constraints::CalculateConstraintOutput<DeviceType>{
            constraints.type,
            constraints.target_node_index,
            constraints.axes,
            state.x0,
            state.q,
            state.v,
            state.vd,
            constraints.output,
        }
    );
    Kokkos::deep_copy(constraints.host_output, constraints.output);

    // If the error was NaN, return false
    return !std::isnan(err);
}

/**
 * @brief Advance solution by one time step with bisection-based load reduction for static analysis
 *
 * For static analysis, uses bisection search to find the maximum load factor that allows
 * convergence. Scales external loads and maintains converged state for rollback on failure.
 * For dynamic analysis, performs a single nonlinear solve.
 *
 * @return true if a converged solution is found at load factor 1.0, otherwise false
 */
template <typename DeviceType>
inline bool Step(
    StepParameters& parameters, Solver<DeviceType>& solver, Elements<DeviceType>& elements,
    State<DeviceType>& state, Constraints<DeviceType>& constraints
) {
    //--------------------------------------------------------------------------
    // Dynamic analysis -> just do a single nonlinear increment
    //--------------------------------------------------------------------------
    if (parameters.is_dynamic_solve) {
        return SolveStep(parameters, solver, elements, state, constraints);
    }

    //--------------------------------------------------------------------------
    // Static analysis -> load reduction strategy to help with convergence
    //--------------------------------------------------------------------------
    const Kokkos::View<double* [6], DeviceType> loads_baseline("loads_baseline", state.f.extent(0));
    Kokkos::deep_copy(loads_baseline, state.f);

    // Load reduction variables initialization
    bool solved{false};
    double load_factor_low{0.};      // lower bound for the load factor
    double load_factor_current{1.};  // current load factor
    double load_factor_high{1.};     // upper bound for the load factor

    //-------------------------------
    // Bisection method
    //-------------------------------
    // Load reduction loop to help with convergence
    auto state_last_converged = CloneState<DeviceType>(state);
    for (size_t j = 0; j < parameters.static_load_retries; ++j) {
        // Scale the external loads
        Kokkos::parallel_for(
            "ScaleLoads",
            Kokkos::RangePolicy<typename DeviceType::execution_space>(0, state.f.extent(0)),
            KOKKOS_LAMBDA(const int i) {
                for (int k = 0; k < 6; ++k) {
                    state.f(i, k) = loads_baseline(i, k) * load_factor_current;
                }
            }
        );

        // TODO Scale the gravity load as well?

        // Attempt a single nonlinear solve at the current load factor
        const bool converged = SolveStep(parameters, solver, elements, state, constraints);

        if (converged) {
            // Record the last converged state - if we're at full load our work is done
            CopyStateData<DeviceType>(state_last_converged, state);
            if (std::abs(load_factor_current - 1.) < 1e-10) {
                solved = true;
                break;
            }
            // We have a successful converged lower bound -> increase lower bound
            // to current load factor + set current load factor to full load
            load_factor_low = load_factor_current;
            load_factor_current = 1.;
        } else {
            // Not converged -> shrink load interval by bisection i.e. upper bounds is reduced to
            // current load factor + new current load factor is set to the average of the bounds i.e.
            // bisected
            load_factor_high = load_factor_current;
            load_factor_current = 0.5 * (load_factor_low + load_factor_high);
            // Roll back to the last converged state
            CopyStateData<DeviceType>(state, state_last_converged);
        }
    }
    return solved;
}

}  // namespace kynema
