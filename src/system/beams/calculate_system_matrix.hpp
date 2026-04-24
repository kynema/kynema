#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculateSystemMatrix {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    ConstView<double* [6][6]> tangent;
    ConstView<size_t**> node_state_indices;
    ConstView<double** [6][6]> stiffness_matrix_terms;
    ConstView<double** [6][6]> inertia_matrix_terms;
    View<double*** [6][6]> system_matrix_terms;

    KOKKOS_FUNCTION
    void operator()(size_t node_12) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using CopyMatrix = KokkosBatched::SerialCopy<>;
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using Gemm = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;

        const auto node_1 = node_12 / num_nodes;
        const auto node_2 = node_12 % num_nodes;
        const auto node_T = node_state_indices(element, node_2);

        auto S_data = Array<double, 36>{};
        auto T_data = Array<double, 36>{};
        auto STpI_data = Array<double, 36>{};

        const auto S = View<double[6][6]>(S_data.data());
        const auto T = View<double[6][6]>(T_data.data());
        const auto STpI = View<double[6][6]>(STpI_data.data());

        CopyMatrix::invoke(subview(stiffness_matrix_terms, node_1, node_2, ALL, ALL), S);
        CopyMatrix::invoke(subview(tangent, node_T, ALL, ALL), T);
        CopyMatrix::invoke(subview(inertia_matrix_terms, node_1, node_2, ALL, ALL), STpI);

        Gemm::invoke(1., S, T, 1., STpI);

        CopyMatrix::invoke(STpI, subview(system_matrix_terms, element, node_1, node_2, ALL, ALL));
    }
};

}  // namespace kynema_fmb::beams
