#pragma once

#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::beams {
template <typename DeviceType>
struct CalculateStrain {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& x0_prime, const ConstView<double[3]>& u_prime,
        const ConstView<double[4]>& r, const ConstView<double[4]>& r_prime,
        const View<double[6]>& strain
    ) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        using Kokkos::Array;

        auto R_x0_prime_data = Array<double, 3>{};
        auto R_x0_prime = View<double[3]>(R_x0_prime_data.data());

        math::RotateVectorByQuaternion(r, x0_prime, R_x0_prime);
        KokkosBlas::serial_axpy(-1., u_prime, R_x0_prime);
        KokkosBlas::serial_axpy(-1., x0_prime, R_x0_prime);

        auto E_data = Array<double, 12>{};
        auto E = View<double[3][4]>(E_data.data());
        math::QuaternionDerivative(r, E);
        auto e2_data = Array<double, 4>{};
        auto e2 = View<double[4]>{e2_data.data()};
        Gemv::invoke(2., E, r_prime, 0., e2);

        strain(0) = -R_x0_prime(0);
        strain(1) = -R_x0_prime(1);
        strain(2) = -R_x0_prime(2);
        strain(3) = e2(0);
        strain(4) = e2(1);
        strain(5) = e2(2);
    }
};
}  // namespace kynema_fmb::beams
