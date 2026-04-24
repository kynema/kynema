#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename shape_matrix_type, typename node_type, typename qp_type>
KOKKOS_INLINE_FUNCTION void InterpVector3(
    const shape_matrix_type& shape_matrix, const node_type& node_v, const qp_type& qp_v
) {
    for (auto qp = 0; qp < qp_v.extent_int(0); ++qp) {
        auto local_total = Kokkos::Array<double, 3>{};
        for (auto node = 0; node < node_v.extent_int(0); ++node) {
            const auto phi = shape_matrix(node, qp);
            for (auto component = 0; component < 3; ++component) {
                local_total[component] += node_v(node, component) * phi;
            }
        }
        for (auto component = 0; component < 3; ++component) {
            qp_v(qp, component) = local_total[component];
        }
    }
}

template <typename shape_matrix_type, typename node_type, typename qp_type>
KOKKOS_INLINE_FUNCTION void InterpVector4(
    const shape_matrix_type& shape_matrix, const node_type& node_v, const qp_type& qp_v
) {
    for (auto qp = 0; qp < qp_v.extent_int(0); ++qp) {
        auto local_total = Kokkos::Array<double, 4>{};
        for (auto node = 0; node < node_v.extent_int(0); ++node) {
            const auto phi = shape_matrix(node, qp);
            for (auto component = 0; component < 4; ++component) {
                local_total[component] += node_v(node, component) * phi;
            }
        }
        for (auto component = 0; component < 4; ++component) {
            qp_v(qp, component) = local_total[component];
        }
    }
}

template <typename shape_matrix_type, typename node_type, typename qp_type>
KOKKOS_INLINE_FUNCTION void InterpQuaternion(
    const shape_matrix_type& shape_matrix, const node_type& node_v, const qp_type& qp_v
) {
    InterpVector4(shape_matrix, node_v, qp_v);
    static constexpr auto length_zero_result = Kokkos::Array<double, 4>{1., 0., 0., 0.};
    for (auto qp = 0; qp < qp_v.extent_int(0); ++qp) {
        auto length = Kokkos::sqrt(
            Kokkos::pow(qp_v(qp, 0), 2) + Kokkos::pow(qp_v(qp, 1), 2) + Kokkos::pow(qp_v(qp, 2), 2) +
            Kokkos::pow(qp_v(qp, 3), 2)
        );
        if (length == 0.) {
            for (auto component = 0; component < 4; ++component) {
                qp_v(qp, component) = length_zero_result[component];
            }
        } else {
            for (auto component = 0; component < 4; ++component) {
                qp_v(qp, component) /= length;
            }
        }
    }
}

template <typename shape_matrix_type, typename jacobian_type, typename node_type, typename qp_type>
KOKKOS_INLINE_FUNCTION void InterpVector3Deriv(
    const shape_matrix_type& shape_matrix_deriv, const jacobian_type& jacobian,
    const node_type& node_v, const qp_type& qp_v
) {
    InterpVector3(shape_matrix_deriv, node_v, qp_v);
    for (auto qp = 0; qp < qp_v.extent_int(0); ++qp) {
        const auto jac = jacobian(qp);
        for (auto component = 0; component < qp_v.extent_int(1); ++component) {
            qp_v(qp, component) /= jac;
        }
    }
}

template <typename shape_matrix_type, typename jacobian_type, typename node_type, typename qp_type>
KOKKOS_INLINE_FUNCTION void InterpVector4Deriv(
    const shape_matrix_type& shape_matrix_deriv, const jacobian_type& jacobian,
    const node_type& node_v, const qp_type& qp_v
) {
    InterpVector4(shape_matrix_deriv, node_v, qp_v);
    for (auto qp = 0; qp < qp_v.extent_int(0); ++qp) {
        const auto jac = jacobian(qp);
        for (auto component = 0; component < qp_v.extent_int(1); ++component) {
            qp_v(qp, component) /= jac;
        }
    }
}

}  // namespace kynema_fmb::beams
