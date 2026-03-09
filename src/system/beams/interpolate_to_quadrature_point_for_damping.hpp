#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema::beams {

template <typename DeviceType>
struct InterpolateToQuadraturePointForDamping {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double jacobian, const ConstLeftView<double*>& shape_interp,
        const ConstLeftView<double*>& shape_deriv, const ConstView<double* [7]>& node_u,
        const ConstView<double* [7]>& node_v, const View<double[3]>& u, const View<double[4]>& r,
        const View<double[3]>& u_prime, const View<double[4]>& r_prime
    ) {
        for (auto node = 0U; node < node_u.extent(0); ++node) {
            const auto phi = shape_interp(node);
            const auto dphiJ = shape_deriv(node) / jacobian;
            for (auto component = 0U; component < 3U; ++component) {
                u(component) += node_u(node, component) * phi;
                u_prime(component) += node_u(node, component) * dphiJ;
            }
            for (auto component = 0U; component < 4U; ++component) {
                r(component) += node_u(node, component + 3) * phi;
                r_prime(component) += node_u(node, component + 3) * dphiJ;
            }
        }
        const auto r_length = KokkosBlas::serial_nrm2(r);
        for (auto component = 0U; component < 4U; ++component) {
            r(component) /= r_length;
        }
    }
};
}  // namespace kynema::beams
