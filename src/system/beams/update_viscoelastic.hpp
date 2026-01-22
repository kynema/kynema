#pragma once

#include <Kokkos_Core.hpp>

namespace kynema::beams {

template <typename DeviceType>
struct UpdateViscoelastic {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double h, double tau_i, const ConstView<double[6]>& strain_dot_n,
        const ConstView<double[6]>& strain_dot_n1, const View<double[6]>& visco_hist
    ) {
        const auto tmp = -1. * h / tau_i;
        const auto tmp_exp = Kokkos::exp(tmp);
        for (auto component = 0; component < 6; ++component) {
            visco_hist(component) = tmp_exp * visco_hist(component) +
                                    (h / 2.) * tmp_exp * strain_dot_n(component) +
                                    (h / 2.) * strain_dot_n1(component);
        }
    }
};

}  // namespace kynema::beams
