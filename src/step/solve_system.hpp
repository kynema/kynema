#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "solver/condition_system.hpp"
#include "solver/linear_solver/dss_numeric.hpp"
#include "solver/linear_solver/dss_solve.hpp"
#include "solver/solver.hpp"
#include "step_parameters.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void SolveSystem(StepParameters& parameters, Solver<DeviceType>& solver) {
    auto region = Kokkos::Profiling::ScopedRegion("Solve System");

    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;

    Kokkos::parallel_for(
        "ConditionR", RangePolicy(0, solver.num_system_dofs),
        solver::ConditionR<DeviceType>{parameters.conditioner, solver.b}
    );

    KokkosBlas::scal(solver.b, -1., solver.b);

    {
        auto solve_region = Kokkos::Profiling::ScopedRegion("Linear Solve");
        dss::numeric_factorization(solver.handle, solver.A);
        dss::solve(solver.handle, solver.A, solver.b, solver.x);
    }

    Kokkos::parallel_for(
        "UnconditionSolution", RangePolicy(0, solver.num_dofs - solver.num_system_dofs),
        solver::UnconditionSolution<DeviceType>{
            solver.num_system_dofs, parameters.conditioner, solver.x
        }
    );
}

}  // namespace kynema_fmb::step
