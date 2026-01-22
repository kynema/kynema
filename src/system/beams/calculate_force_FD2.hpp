#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateForceFD2 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[6][6]>& D, const ConstView<double[3]>& xr_prime, const ConstView<double[3]>& u_prime, const ConstView<double[6]>& strain_dot, const View<double[6]>& FD2) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Transpose = KokkosBlas::Trans::Transpose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        using GemvT = KokkosBlas::SerialGemv<Transpose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto xr_plus_u_prime_data = Kokkos::Array<double, 3>{};
        auto xr_plus_u_prime = View<double[3]>(xr_plus_u_prime_data.data());
        for (auto component = 0; component < 3; ++component) {
            xr_plus_u_prime(component) = u_prime(component) + xr_prime(component);
        }

        auto xr_plus_u_prime_tilde_data = Kokkos::Array<double, 9>{};
        auto xr_plus_u_prime_tilde = View<double[3][3]>(xr_plus_u_prime_tilde_data.data());

        math::VecTilde(xr_plus_u_prime, xr_plus_u_prime_tilde);

        for (auto component = 0; component < 3; ++component) {
            FD2(component) = 0.;
        }

        auto eps_kappa_sum_data = Kokkos::Array<double, 3>{};
        auto eps_kappa_sum = View<double[3]>(eps_kappa_sum_data.data());

        auto eps = subview(D, make_pair(0, 3), make_pair(0, 3));
        Gemv::invoke(1., eps, subview(strain_dot, make_pair(0, 3)), 0., eps_kappa_sum);
        auto kappa = subview(D, make_pair(0, 3), make_pair(3, 6));
        Gemv::invoke(1., kappa, subview(strain_dot, make_pair(3, 6)), 1., eps_kappa_sum);
        GemvT::invoke(1., xr_plus_u_prime_tilde, eps_kappa_sum, 0., subview(FD2, make_pair(3, 6)));
    }
};

}
