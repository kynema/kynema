#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which applies the given factor to the system RHS vector
 */
template <typename DeviceType>
struct ConditionR {
    double conditioner;
    Kokkos::View<double* [1], Kokkos::LayoutLeft, DeviceType> R;

    KOKKOS_FUNCTION
    void operator()(int i) const { R(i, 0) *= conditioner; }
};

/**
 * @brief A Kernel which divides the RHS vector terms corresponding to the constraints
 * by a given conditioner factor
 */
template <typename DeviceType>
struct UnconditionSolution {
    size_t num_system_dofs;
    double conditioner;
    Kokkos::View<double* [1], Kokkos::LayoutLeft, DeviceType> x;

    KOKKOS_FUNCTION
    void operator()(size_t i) const { x(i + num_system_dofs, 0) /= conditioner; }
};

}  // namespace kynema_fmb::solver
