#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema_fmb::system {

/**
 * @brief A Kernel for calculating the tangent operator at each node
 */
template <typename DeviceType>
struct CalculateTangentOperator {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    double h;
    ConstView<double* [6]> q_delta;
    View<double* [6][6]> T_gbl;

    KOKKOS_FUNCTION
    void operator()(int node) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;
        using CopyMatrix = KokkosBatched::SerialCopy<>;
        using CopyMatrixTranspose = KokkosBatched::SerialCopy<KokkosBatched::Trans::Transpose>;
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using Gemm = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;

        auto T_data = Array<double, 36>{};
        const auto T = View<double[6][6]>(T_data.data());

        auto rv_data = Array<double, 3>{};
        auto rv = View<double[3]>{rv_data.data()};

        auto m1_data = Array<double, 9>{};
        auto m1 = View<double[3][3]>(m1_data.data());
        auto m2_data = Array<double, 9>{};
        auto m2 = View<double[3][3]>(m2_data.data());

        for (auto k = 0; k < 6; ++k) {
            T(k, k) = 1.;
        }

        KokkosBlas::serial_axpy(h, subview(q_delta, node, make_pair(3, 6)), rv);
        auto phi = KokkosBlas::serial_nrm2(rv);
        const auto tmp1 = (phi > 1.e-16) ? (Kokkos::cos(phi) - 1.) / (phi * phi) : 0.;
        const auto tmp2 = (phi > 1.e-16) ? (1. - Kokkos::sin(phi) / phi) / (phi * phi) : 0.;

        math::VecTilde(rv, m2);
        CopyMatrix::invoke(m2, m1);
        Gemm::invoke(tmp2, m2, m2, tmp1, m1);

        KokkosBlas::serial_axpy(1., m1, subview(T, make_pair(3, 6), make_pair(3, 6)));

        CopyMatrixTranspose::invoke(T, subview(T_gbl, node, ALL, ALL));
    }
};

}  // namespace kynema_fmb::system
