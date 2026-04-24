#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::state {

/**
 * @brief A Kernel to predict the next state values before nonlinear iteration begins for
 * a given time step
 */
template <typename DeviceType>
struct CalculateNextState {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;

    double h;
    double alpha_f;
    double alpha_m;
    double beta;
    double gamma;
    View<double* [6]> q_delta;
    View<double* [6]> v;
    View<double* [6]> vd;
    View<double* [6]> a;

    KOKKOS_FUNCTION
    void operator()(size_t node) const {
        for (auto component = 0U; component < 6U; ++component) {
            const double v_p = v(node, component);
            const double vd_p = vd(node, component);
            const double a_p = a(node, component);
            vd(node, component) = 0.;
            a(node, component) = (alpha_f * vd_p - alpha_m * a_p) / (1. - alpha_m);
            v(node, component) = v_p + h * (1. - gamma) * a_p + gamma * h * a(node, component);
            q_delta(node, component) = v_p + (0.5 - beta) * h * a_p + beta * h * a(node, component);
        }
    }
};

}  // namespace kynema_fmb::state
