#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "elements/elements.hpp"
#include "solver/contribute_beams_to_sparse_matrix.hpp"
#include "solver/contribute_masses_to_sparse_matrix.hpp"
#include "solver/contribute_springs_to_sparse_matrix.hpp"
#include "solver/solver.hpp"
#include "step_parameters.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void AssembleSystemMatrix(
    StepParameters& parameters, Solver<DeviceType>& solver, Elements<DeviceType>& elements
) {
    auto region = Kokkos::Profiling::ScopedRegion("Assemble System Matrix");

    using TeamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>;

    const auto num_nodes = static_cast<int>(elements.beams.max_elem_nodes);
    const auto num_beams = static_cast<int>(elements.beams.num_elems);
    const auto num_masses = static_cast<int>(elements.masses.num_elems);
    const auto num_springs = static_cast<int>(elements.springs.num_elems);

    const auto vector_length = std::min(num_nodes, TeamPolicy::vector_length_max());
    auto beams_sparse_matrix_policy = TeamPolicy(num_beams, Kokkos::AUTO(), vector_length);
    auto masses_sparse_matrix_policy = TeamPolicy(num_masses, Kokkos::AUTO());
    auto springs_sparse_matrix_policy = TeamPolicy(num_springs, Kokkos::AUTO());

    Kokkos::parallel_for(
        "ContributeBeamsToSparseMatrix", beams_sparse_matrix_policy,
        solver::ContributeBeamsToSparseMatrix<typename Solver<DeviceType>::CrsMatrixType>{
            parameters.conditioner, elements.beams.num_nodes_per_element,
            elements.beams.element_freedom_signature, elements.beams.element_freedom_table,
            elements.beams.system_matrix_terms, solver.A
        }
    );
    Kokkos::parallel_for(
        "ContributeMassesToSparseMatrix", masses_sparse_matrix_policy,
        solver::ContributeMassesToSparseMatrix<typename Solver<DeviceType>::CrsMatrixType>{
            parameters.conditioner, elements.masses.element_freedom_signature,
            elements.masses.element_freedom_table, elements.masses.system_matrix_terms, solver.A
        }
    );
    Kokkos::parallel_for(
        "ContributeSpringsToSparseMatrix", springs_sparse_matrix_policy,
        solver::ContributeSpringsToSparseMatrix<typename Solver<DeviceType>::CrsMatrixType>{
            parameters.conditioner, elements.springs.element_freedom_signature,
            elements.springs.element_freedom_table, elements.springs.stiffness_matrix_terms, solver.A
        }
    );
}

}  // namespace kynema_fmb::step
