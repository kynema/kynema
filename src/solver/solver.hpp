#pragma once

#include <vector>

#include <KokkosSparse.hpp>
#include <Kokkos_Core.hpp>

#include "compute_num_system_dofs.hpp"
#include "create_full_matrix.hpp"
#include "linear_solver/dss_handle.hpp"
#include "linear_solver/dss_symbolic.hpp"

#include "KYNEMA_FMB_config.h"

namespace kynema_fmb {

/** @brief This object manages the assembly and solution of linear system arising from the
 * generalized-alpha based time integration of the dynamic structural problem.
 */
template <typename DeviceType>
struct Solver {
    // Define some types for the solver to make the code more readable
    using ValueType = double;
#ifdef KOKKOS_ENABLE_CUDA
    static constexpr bool use_device =
        std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value;
#if defined(KYNEMA_FMB_ENABLE_CUDSS)
    static constexpr dss::Algorithm algorithm_device = dss::Algorithm::CUDSS;
#elif defined(KYNEMA_FMB_ENABLE_CUSOLVERSP)
    static constexpr dss::Algorithm algorithm_device = dss::Algorithm::CUSOLVER_SP;
#else
    static constexpr dss::Algorithm algorithm_device = dss::Algorithm::NONE;
#endif
#else
    static constexpr bool use_device = false;
    static constexpr dss::Algorithm algorithm_device = dss::Algorithm::NONE;
#endif

#if defined(KYNEMA_FMB_ENABLE_KLU)
    static constexpr dss::Algorithm algorithm_host = dss::Algorithm::KLU;
#elif defined(KYNEMA_FMB_ENABLE_SUPERLU)
    static constexpr dss::Algorithm algorithm_host = dss::Algorithm::SUPERLU;
#elif defined(KYNEMA_FMB_ENABLE_MKL)
    static constexpr dss::Algorithm algorithm_host = dss::Algorithm::MKL;
#elif defined(KYNEMA_FMB_ENABLE_SUPERLU_MT)
    static constexpr dss::Algorithm algorithm_host = dss::Algorithm::SUPERLU_MT;
#elif defined(KYNEMA_FMB_ENABLE_UMFPACK)
    static constexpr dss::Algorithm algorithm_host = dss::Algorithm::UMFPACK;
#else
    static constexpr dss::Algorithm algorithm_host = dss::Algorithm::NONE;
#endif

    static constexpr dss::Algorithm algorithm =
        (use_device && algorithm_device != dss::Algorithm::NONE) ? algorithm_device : algorithm_host;

    static_assert(algorithm != dss::Algorithm::NONE);

    using HandleType = dss::Handle<algorithm>;
#if defined(KYNEMA_FMB_ENABLE_MKL)
    using IndexType = std::conditional<algorithm == dss::Algorithm::MKL, MKL_INT, int>::type;
#else
    using IndexType = int;
#endif

    using CrsMatrixType = KokkosSparse::CrsMatrix<ValueType, IndexType, DeviceType, void, IndexType>;
    using MultiVectorType = Kokkos::View<ValueType* [1], Kokkos::LayoutLeft, DeviceType>;

    template <typename value_type>
    using View = Kokkos::View<value_type, DeviceType>;
    template <typename value_type>
    using ConstView = typename View<value_type>::const_type;

    size_t num_system_nodes;  //< Number of system nodes
    size_t num_system_dofs;   //< Number of system degrees of freedom
    size_t num_dofs;          //< Number of degrees of freedom

    CrsMatrixType A;    //< System matrix
    MultiVectorType b;  //< System RHS
    MultiVectorType x;  //< System solution

    HandleType handle;  //< Handle for internal information needed for the selected linear solver

    std::vector<double> convergence_err;  //< The convergence history of the solver

    /** @brief Constructs the sparse matrix structure for the provided connectivity information
     * and performs the sparse direct solver's symbolic factorization step, which initializes
     * its internal data structures.
     *
     * @param node_IDs View containing the global IDs for each node
     * @param active_dofs View containing the number of active DoFs for each node
     * @param node_freedom_map_table View mapping node indices to DOF indices
     * @param num_nodes_per_element View containing number of nodes per element
     * @param node_state_indices View containing element-to-node connectivity
     * @param num_constraint_dofs Number of constraint degrees of freedom
     * @param base_active_dofs View containing the number of DoFs at the base node
     * @param target_active_dofs View containing the number of DoFs at the target node
     * @param constraint_base_node_freedom_table View containing base node DOFs for constraints
     * @param constraint_target_node_freedom_table View containing target node DOFs for constraints
     * @param constraint_row_range View containing row ranges for each constraint
     */
    Solver(
        const ConstView<size_t*>& node_IDs, const ConstView<size_t*>& active_dofs,
        const ConstView<size_t*>& node_freedom_map_table,
        const ConstView<size_t*>& num_nodes_per_element,
        const ConstView<size_t**>& node_state_indices, size_t num_constraint_dofs,
        const ConstView<size_t*>& base_active_dofs, const ConstView<size_t*>& target_active_dofs,
        const ConstView<size_t* [6]>& constraint_base_node_freedom_table,
        const ConstView<size_t* [6]>& constraint_target_node_freedom_table,
        const ConstView<Kokkos::pair<size_t, size_t>*>& constraint_row_range
    )
        : num_system_nodes(node_IDs.extent(0)),
          num_system_dofs(solver::ComputeNumSystemDofs<DeviceType>(active_dofs)),
          num_dofs(num_system_dofs + num_constraint_dofs),
          A(solver::CreateFullMatrix<CrsMatrixType>::invoke(
              num_system_dofs, num_dofs, base_active_dofs, target_active_dofs,
              constraint_base_node_freedom_table, constraint_target_node_freedom_table,
              constraint_row_range, active_dofs, node_freedom_map_table, num_nodes_per_element,
              node_state_indices
          )),
          b(Kokkos::view_alloc("b", Kokkos::WithoutInitializing), num_dofs),
          x(Kokkos::view_alloc("x", Kokkos::WithoutInitializing), num_dofs) {
        dss::symbolic_factorization(handle, A);
    }
};

}  // namespace kynema_fmb
