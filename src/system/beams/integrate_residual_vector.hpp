#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct IntegrateResidualVectorElement {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    size_t element;
    size_t num_qps;
    ConstView<double*> qp_weight_;
    ConstView<double*> qp_jacobian_;
    ConstLeftView<double**> shape_interp_;
    ConstLeftView<double**> shape_deriv_;
    ConstView<double* [6]> node_FX_;
    ConstView<double* [6]> qp_FE1_;
    ConstView<double* [6]> qp_FE2_;
    ConstView<double* [6]> qp_FD1_;
    ConstView<double* [6]> qp_FD2_;
    ConstView<double* [6]> qp_Fi_;
    ConstView<double* [6]> qp_Fe_;
    ConstView<double* [6]> qp_Fg_;
    View<double** [6]> residual_vector_terms_;

    KOKKOS_FUNCTION
    void operator()(size_t node) const {
        auto local_residual = Kokkos::Array<double, 6>{};
        for (auto qp = 0U; qp < num_qps; ++qp) {
            const auto weight = qp_weight_(qp);
            const auto coeff_c = weight * shape_deriv_(node, qp);
            const auto coeff_dig = weight * qp_jacobian_(qp) * shape_interp_(node, qp);
            for (auto component = 0U; component < 6U; ++component) {
                local_residual[component] +=
                    coeff_c * (qp_FE1_(qp, component) + qp_FD1_(qp, component)) +
                    coeff_dig *
                        (qp_FE2_(qp, component) + qp_FD2_(qp, component) + qp_Fi_(qp, component) -
                         qp_Fe_(qp, component) - qp_Fg_(qp, component));
            }
        }
        for (auto component = 0U; component < 6U; ++component) {
            residual_vector_terms_(element, node, component) =
                local_residual[component] - node_FX_(node, component);
        }
    }
};

}  // namespace kynema_fmb::beams
