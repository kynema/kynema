#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

namespace kynema::beams {

template <typename DeviceType>
struct IntegrateStiffnessMatrixElement {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    size_t num_qps;
    ConstView<double*> qp_weight_;
    ConstView<double*> qp_jacobian_;
    ConstLeftView<double**> shape_interp_;
    ConstLeftView<double**> shape_deriv_;
    ConstView<double* [6][6]> qp_Kuu_;
    ConstView<double* [6][6]> qp_Puu_;
    ConstView<double* [6][6]> qp_Cuu_;
    ConstView<double* [6][6]> qp_Ouu_;
    ConstView<double* [6][6]> qp_Quu_;
    View<double** [6][6]> gbl_M_;

    KOKKOS_FUNCTION
    void operator()(size_t node_simd_node) const {
        using simd_type = Kokkos::Experimental::simd<double>;
        using tag_type = Kokkos::Experimental::vector_aligned_tag;
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        constexpr auto width = simd_type::size();
        const auto extra_component = num_nodes % width == 0U ? 0U : 1U;
        const auto simd_nodes = (num_nodes / width) + extra_component;
        const auto node = node_simd_node / simd_nodes;
        const auto simd_node = (node_simd_node % simd_nodes) * width;

        auto local_M = Array<simd_type, 36>{};

        const auto qp_Kuu = ConstView<double* [36]>(qp_Kuu_.data(), num_qps);
        const auto qp_Puu = ConstView<double* [36]>(qp_Puu_.data(), num_qps);
        const auto qp_Cuu = ConstView<double* [36]>(qp_Cuu_.data(), num_qps);
        const auto qp_Ouu = ConstView<double* [36]>(qp_Ouu_.data(), num_qps);
        const auto qp_Quu = ConstView<double* [36]>(qp_Quu_.data(), num_qps);

        for (auto qp = 0U; qp < num_qps; ++qp) {
            const auto w = simd_type(qp_weight_(qp));
            const auto jacobian = simd_type(qp_jacobian_(qp));
            const auto phi_1 = simd_type(shape_interp_(node, qp));
            auto phi_2 = simd_type{};
            phi_2.copy_from(&shape_interp_(simd_node, qp), tag_type());
            const auto phi_prime_1 = simd_type(shape_deriv_(node, qp));
            auto phi_prime_2 = simd_type{};
            phi_prime_2.copy_from(&shape_deriv_(simd_node, qp), tag_type());
            const auto K = (phi_1 * phi_2) * (w * jacobian);
            const auto P = (phi_1 * phi_prime_2) * w;
            const auto C = (phi_prime_1 * phi_prime_2) * (w / jacobian);
            const auto O = (phi_prime_1 * phi_2) * w;
            const auto Kuu_local = subview(qp_Kuu, qp, ALL);
            const auto Quu_local = subview(qp_Quu, qp, ALL);
            const auto Puu_local = subview(qp_Puu, qp, ALL);
            const auto Cuu_local = subview(qp_Cuu, qp, ALL);
            const auto Ouu_local = subview(qp_Ouu, qp, ALL);
            for (auto component = 0; component < 36; ++component) {
                const auto Kuu = simd_type(Kuu_local(component));
                const auto Quu = simd_type(Quu_local(component));
                const auto Puu = simd_type(Puu_local(component));
                const auto Cuu = simd_type(Cuu_local(component));
                const auto Ouu = simd_type(Ouu_local(component));
                local_M[component] =
                    local_M[component] + K * (Kuu + Quu) + P * Puu + C * Cuu + O * Ouu;
            }
        }

        const auto num_lanes = Kokkos::min(width, num_nodes - simd_node);
        const auto global_M = View<double** [36]>(gbl_M_.data(), num_nodes, num_nodes);
        const auto M_slice =
            subview(global_M, node, make_pair(simd_node, simd_node + num_lanes), ALL);

        for (auto lane = 0U; lane < num_lanes; ++lane) {
            for (auto component = 0; component < 36; ++component) {
                M_slice(lane, component) = local_M[component][lane];
            }
        }
    }
};

}  // namespace kynema::beams
