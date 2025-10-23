#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

namespace kynema::beams {

template <typename DeviceType>
struct IntegrateInertiaMatrixElement {
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
    ConstView<double* [6][6]> qp_Muu_;
    ConstView<double* [6][6]> qp_Guu_;
    double beta_prime_;
    double gamma_prime_;
    Kokkos::View<double** [6][6], DeviceType> gbl_M_;

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
        const auto num_simd_nodes = (num_nodes / width) + extra_component;
        const auto node = node_simd_node / num_simd_nodes;
        const auto simd_node = (node_simd_node % num_simd_nodes) * width;

        auto local_M = Array<simd_type, 36>{};

        const auto qp_Muu = ConstView<double* [36]>(qp_Muu_.data(), num_qps);
        const auto qp_Guu = ConstView<double* [36]>(qp_Guu_.data(), num_qps);

        for (auto qp = 0U; qp < num_qps; ++qp) {
            const auto w = simd_type(qp_weight_(qp));
            const auto jacobian = simd_type(qp_jacobian_(qp));
            const auto phi_1 = simd_type(shape_interp_(node, qp));
            auto phi_2 = simd_type{};
            phi_2.copy_from(&shape_interp_(simd_node, qp), tag_type());
            const auto coeff = phi_1 * phi_2 * w * jacobian;
            const auto Muu_local = subview(qp_Muu, qp, ALL);
            const auto Guu_local = subview(qp_Guu, qp, ALL);
            for (auto component = 0; component < 36; ++component) {
                const auto contribution = simd_type(
                    (beta_prime_ * Muu_local(component)) + (gamma_prime_ * Guu_local(component))
                );
                local_M[component] = local_M[component] + coeff * contribution;
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
