#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculateTemporaryVariables {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& x0_prime, const ConstView<double[3]>& u_prime,
        const View<double[3][3]>& x0pupSS
    ) {
        using Kokkos::Array;

        auto x0pup_data = Array<double, 3>{x0_prime(0), x0_prime(1), x0_prime(2)};
        const auto x0pup = View<double[3]>(x0pup_data.data());

        KokkosBlas::serial_axpy(1., u_prime, x0pup);
        math::VecTilde(x0pup, x0pupSS);
    }
};
}  // namespace kynema_fmb::beams
