#pragma once

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

namespace kynema_fmb::beams::tests {

template <typename ValueType, typename DataType>
typename Kokkos::View<ValueType>::const_type CreateView(
    const std::string& name, const DataType& data
) {
    const auto view = Kokkos::View<ValueType>(Kokkos::view_alloc(name, Kokkos::WithoutInitializing));
    const auto host = typename Kokkos::View<ValueType, Kokkos::HostSpace>::const_type(data.data());
    const auto mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, view);
    Kokkos::deep_copy(mirror, host);
    Kokkos::deep_copy(view, mirror);
    return view;
}

inline auto create_shape_interp_OneNodeOneQP() {
    return CreateView<double[1][1][1]>("shape_interp", std::array<double, 1>{2.});
}

inline auto create_shape_deriv_OneNodeOneQP() {
    return CreateView<double[1][1][1]>("shape_deriv", std::array<double, 1>{4.});
}

inline auto create_jacobian_OneQP() {
    return CreateView<double[1][1]>("jacobian", std::array<double, 1>{2.});
}

inline auto create_shape_interp_OneNodeTwoQP() {
    return CreateView<double[1][1][2]>("shape_interp", std::array{2., 4.});
}

inline auto create_shape_deriv_OneNodeTwoQP() {
    return CreateView<double[1][1][2]>("shape_deriv", std::array{4., 9.});
}

inline auto create_jacobian_TwoQP() {
    return CreateView<double[1][2]>("jacobian", std::array{2., 3.});
}

inline auto create_shape_interp_TwoNodeTwoQP() {
    return CreateView<double[1][2][2]>("shape_interp", std::array{2., 3., 4., 5.});
}

inline auto create_shape_deriv_TwoNodeTwoQP() {
    return CreateView<double[1][2][2]>("shape_deriv", std::array{4., 9., 8., 15.});
}

}  // namespace kynema_fmb::beams::tests
