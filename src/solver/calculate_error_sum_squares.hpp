#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief Reduction Kernel which calculates the sum of the square of the errors for each
 * node in the system for use in computing the system convergence.
 */
template <typename DeviceType>
struct CalculateSystemErrorSumSquares {
    using value_type = double;
    using TeamPolicy = typename Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    double atol;
    double rtol;
    double h;
    ConstView<size_t*> active_dofs;
    ConstView<size_t*> node_freedom_map_table;
    ConstView<double* [6]> q_delta;
    ConstLeftView<double* [1]> x;

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t node, double& err) const {
        const auto n_node_dofs = active_dofs(node);
        const auto node_first_dof = node_freedom_map_table(node);
        for (auto component = 0U; component < n_node_dofs; ++component) {
            const auto pi = x(node_first_dof + component, 0);
            const auto xi = q_delta(node, component) * h;
            const auto err_sqrt = pi / (atol + rtol * Kokkos::abs(xi));
            err += err_sqrt * err_sqrt;
        }
    }
};

/**
 * @brief Reduction Kernel which calculates the sum of the squares of the error
 * for each constraint for use in computing system convergence
 */
template <typename DeviceType>
struct CalculateConstraintsErrorSumSquares {
    using value_type = double;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    double atol;
    double rtol;
    size_t num_system_dofs;
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    ConstView<double* [6]> lambda;
    ConstLeftView<double* [1]> x;

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t constraint, double& err) const {
        const auto first_index = row_range(constraint).first;
        const auto max_index = row_range(constraint).second;
        for (auto row = first_index; row < max_index; ++row) {
            const auto pi = x(num_system_dofs + row, 0);
            const auto xi = lambda(constraint, row - first_index);
            const auto err_sqrt = pi / (atol + rtol * Kokkos::abs(xi));
            err += err_sqrt * err_sqrt;
        }
    }
};

}  // namespace kynema_fmb::solver
