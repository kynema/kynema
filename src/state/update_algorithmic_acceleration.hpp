#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::state {

/**
 * @brief A Kernel to update the algorithmic acceleration based on the acceleration and
 * generalized alpha solver parameters
 */
template <typename DeviceType>
struct UpdateAlgorithmicAcceleration {
    typename Kokkos::View<double* [6], DeviceType> acceleration;
    typename Kokkos::View<double* [6], DeviceType>::const_type vd;
    double alpha_f;
    double alpha_m;

    KOKKOS_FUNCTION
    void operator()(int node) const {
        for (auto component = 0; component < 6; ++component) {
            acceleration(node, component) += (1. - alpha_f) / (1. - alpha_m) * vd(node, component);
        }
    }
};

}  // namespace kynema_fmb::state
