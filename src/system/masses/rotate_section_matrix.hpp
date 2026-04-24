#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::masses {

template <typename DeviceType>
struct RotateSectionMatrix {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[4]>& xr, const ConstView<double[6][6]>& Cstar,
        const View<double[6][6]>& Cuu
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using GemmNT = KokkosBatched::SerialGemm<NoTranspose, Transpose, Default>;
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto RR0_data = Array<double, 9>{};
        auto RR0 = View<double[3][3]>(RR0_data.data());
        math::QuaternionToRotationMatrix(xr, RR0);

        const auto Cstar_top = subview(Cstar, make_pair(0, 3), ALL);
        const auto Cstar_bottom = subview(Cstar, make_pair(3, 6), ALL);
        const auto Cuu_left = subview(Cuu, ALL, make_pair(0, 3));
        const auto Cuu_right = subview(Cuu, ALL, make_pair(3, 6));

        auto ctmp_data = Array<double, 36>{};
        const auto Ctmp = View<double[6][6]>(ctmp_data.data());
        const auto Ctmp_top = subview(Ctmp, make_pair(0, 3), ALL);
        const auto Ctmp_bottom = subview(Ctmp, make_pair(3, 6), ALL);
        const auto Ctmp_left = subview(Ctmp, ALL, make_pair(0, 3));
        const auto Ctmp_right = subview(Ctmp, ALL, make_pair(3, 6));

        GemmNN::invoke(1., RR0, Cstar_top, 0., Ctmp_top);
        GemmNN::invoke(1., RR0, Cstar_bottom, 0., Ctmp_bottom);
        GemmNT::invoke(1., Ctmp_left, RR0, 0., Cuu_left);
        GemmNT::invoke(1., Ctmp_right, RR0, 0., Cuu_right);
    }
};
}  // namespace kynema_fmb::masses
