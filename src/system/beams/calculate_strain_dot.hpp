#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"

namespace kynema::beams {
template <typename DeviceType>
struct CalculateStrainDot {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[3]>& kappa, const ConstView<double[4]>& r, const ConstView<double[3]>& omega, const ConstView<double[3]>& u_dot_prime, const ConstView<double[3]>& omega_prime, const ConstView<double[3]>& xr_prime, const View<double[6]>& strain_dot) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Transpose = KokkosBlas::Trans::Transpose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        using GemvT = KokkosBlas::SerialGemv<Transpose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto omega_tilde_data = Kokkos::Array<double, 9>{};
        auto omega_tilde = View<double[3][3]>(omega_tilde_data.data());
        math::VecTilde(omega, omega_tilde);

        auto r_xr_prime_data = Kokkos::Array<double, 3>{};
        auto r_xr_prime = View<double[3]>(r_xr_prime_data.data());
        math::RotateVectorByQuaternion(r, xr_prime, r_xr_prime);

        for (auto component = 0; component < 3; ++component) {
            strain_dot(component) = u_dot_prime(component);
        }
        Gemv(1., omega_tilde, r_xr_prime, 1., subview(strain_dot, make_pair(0, 3)));
        for (auto component = 0; component < 3; ++component) {
            strain_dot(component + 3) = omega_prime(component);
        }
        Gemv(1., omega_tilde, kappa, 1., subview(strain_dot, make_pair(3, 6)));
    }
};

}
