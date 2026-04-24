#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct InterpolateToQuadraturePointForInertia {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstLeftView<double*>& shape_interp, const ConstView<double* [7]>& node_u,
        const ConstView<double* [6]>& node_u_dot, const ConstView<double* [6]>& node_u_ddot,
        const View<double[4]>& r, const View<double[3]>& u_ddot, const View<double[3]>& omega,
        const View<double[3]>& omega_dot
    ) {
        for (auto node = 0U; node < node_u.extent(0); ++node) {
            const auto phi = shape_interp(node);
            for (auto component = 0U; component < 3U; ++component) {
                u_ddot(component) += node_u_ddot(node, component) * phi;
                omega(component) += node_u_dot(node, component + 3U) * phi;
                omega_dot(component) += node_u_ddot(node, component + 3U) * phi;
            }
            for (auto component = 0U; component < 4U; ++component) {
                r(component) += node_u(node, component + 3) * phi;
            }
        }
        const auto r_length = KokkosBlas::serial_nrm2(r);
        for (auto component = 0U; component < 4U; ++component) {
            r(component) /= r_length;
        }
    }
};
}  // namespace kynema_fmb::beams
