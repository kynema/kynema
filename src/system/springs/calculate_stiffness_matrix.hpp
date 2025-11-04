#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema::springs {

template <typename DeviceType>
struct CalculateStiffnessMatrix {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double c1, double c2, const ConstView<double[3]>& r, double l, const View<double[3][3]>& a
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using GemmDefault = KokkosBatched::Algo::Gemm::Default;
        using Gemm = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, GemmDefault>;

        auto r_tilde_data = Kokkos::Array<double, 9>{};
        auto r_tilde = View<double[3][3]>(r_tilde_data.data());

        const double diag_term = c1 - (c2 * l * l);
        a(0, 0) = diag_term;
        a(1, 1) = diag_term;
        a(2, 2) = diag_term;

        math::VecTilde(r, r_tilde);
        Gemm::invoke(-c2, r_tilde, r_tilde, 1., a);
    }
};
}  // namespace kynema::springs
