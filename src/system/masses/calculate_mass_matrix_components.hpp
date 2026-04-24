#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::masses {

template <typename DeviceType>
KOKKOS_INLINE_FUNCTION void CalculateEta(
    const typename Kokkos::View<double[6][6], DeviceType>::const_type& Muu,
    const Kokkos::View<double[3], DeviceType>& eta
) {
    const auto m = Muu(0, 0);

    eta(0) = Muu(5, 1) / m;
    eta(1) = -Muu(5, 0) / m;
    eta(2) = Muu(4, 0) / m;
}

template <typename DeviceType>
KOKKOS_INLINE_FUNCTION void CalculateRho(
    const typename Kokkos::View<double[6][6], DeviceType>::const_type& Muu,
    const Kokkos::View<double[3][3], DeviceType>& rho
) {
    KokkosBatched::SerialCopy<>::invoke(
        Kokkos::subview(Muu, Kokkos::make_pair(3, 6), Kokkos::make_pair(3, 6)), rho
    );
}

}  // namespace kynema_fmb::masses
