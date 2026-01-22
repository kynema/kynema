#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema::beams {

template <typename DeviceType>
struct CalculateForceFD2Viscoelastic {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double h, double tau_i, const ConstView<double[6][6]>& kv_i,
        const ConstView<double[6][6]>& rr0, const ConstView<double[6]>& strain_dot_n,
        const ConstView<double[6]>& strain_dot_n1, const ConstView<double[6]>& visco_hist,
        const View<double[6]>& FD1
    ) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Transpose = KokkosBlas::Trans::Transpose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        using GemvT = KokkosBlas::SerialGemv<Transpose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        const auto tmp = -1. * h / tau_i;
        const auto tmp_exp = Kokkos::exp(tmp);
        auto visco_curr_data = Kokkos::Array<double, 6>{};
        auto visco_curr = View<double[6]>(visco_curr_data.data());

        for (auto component = 0; component < 6; ++component) {
            visco_curr(component) = tmp_exp * visco_hist(component) +
                                    (h / 2.) * tmp_exp * strain_dot_n(component) +
                                    (h / 2.) * strain_dot_n1(component);
        }

        auto fd_tmp_data = Kokkos::Array<double, 6>{};
        auto fd_tmp = View<double[6]>(fd_tmp_data.data());

        Gemv(1., kv_i, visco_curr, 0., fd_tmp);

        Gemv(1., rr0, fd_tmp, 0., FD1);
    }
};

}  // namespace kynema::beams
