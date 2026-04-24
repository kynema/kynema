#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::constraints {

/**
 * @brief A Kernel to update the predicted constraint Lagrange multiplier values at each
 * nonlinear iteration
 */
template <typename DeviceType>
struct UpdateLambdaPrediction {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    size_t num_system_dofs;
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    ConstLeftView<double* [1]> x;
    View<double* [6]> lambda;

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        const auto first_index = row_range(constraint).first;
        const auto max_index = row_range(constraint).second;
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        for (auto row = first_index; row < max_index; ++row) {
            if constexpr (force_atomic) {
                Kokkos::atomic_add(
                    &lambda(constraint, row - first_index), x(num_system_dofs + row, 0)
                );
            } else {
                lambda(constraint, row - first_index) += x(num_system_dofs + row, 0);
            }
        }
    }
};

}  // namespace kynema_fmb::constraints
