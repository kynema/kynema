#pragma once

#include <Kokkos_Core.hpp>

#include "system/springs/calculate_distance_components.hpp"
#include "system/springs/calculate_force_coefficients.hpp"
#include "system/springs/calculate_force_vectors.hpp"
#include "system/springs/calculate_length.hpp"
#include "system/springs/calculate_stiffness_matrix.hpp"

namespace kynema_fmb::springs {

template <typename DeviceType>
struct CalculateQuadraturePointValues {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<double* [7]> Q;

    ConstView<size_t* [2]> node_state_indices;
    ConstView<double* [3]> x0_;
    ConstView<double*> l_ref_;
    ConstView<double*> k_;

    View<double* [2][3]> residual_vector_terms;
    View<double* [2][2][3][3]> stiffness_matrix_terms;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        using Kokkos::Array;

        const auto index_0 = node_state_indices(element, 0);
        const auto index_1 = node_state_indices(element, 1);

        const auto x0_data = Array<double, 3>{x0_(element, 0), x0_(element, 1), x0_(element, 2)};
        const auto u1_data = Array<double, 3>{Q(index_0, 0), Q(index_0, 1), Q(index_0, 2)};
        const auto u2_data = Array<double, 3>{Q(index_1, 0), Q(index_1, 1), Q(index_1, 2)};
        auto r_data = Array<double, 3>{};
        auto f_data = Array<double, 3>{};
        auto a_data = Array<double, 9>{};

        const auto x0 = ConstView<double[3]>(x0_data.data());
        const auto u1 = ConstView<double[3]>(u1_data.data());
        const auto u2 = ConstView<double[3]>(u2_data.data());
        const auto r = View<double[3]>(r_data.data());
        const auto f = View<double[3]>(f_data.data());
        const auto a = View<double[3][3]>(a_data.data());

        const auto l_ref = l_ref_(element);
        const auto k = k_(element);

        springs::CalculateDistanceComponents<DeviceType>::invoke(x0, u1, u2, r);
        const auto l = springs::CalculateLength<DeviceType>(r);
        const auto c1 = springs::CalculateForceCoefficient1<DeviceType>(k, l_ref, l);
        const auto c2 = springs::CalculateForceCoefficient2<DeviceType>(k, l_ref, l);
        springs::CalculateForceVectors<DeviceType>::invoke(r, c1, f);
        springs::CalculateStiffnessMatrix<DeviceType>::invoke(c1, c2, r, l, a);

        for (auto component = 0; component < 3; ++component) {
            residual_vector_terms(element, 0, component) = f(component);
            residual_vector_terms(element, 1, component) = -f(component);
        }

        for (auto component_1 = 0; component_1 < 3; ++component_1) {
            for (auto component_2 = 0; component_2 < 3; ++component_2) {
                stiffness_matrix_terms(element, 0, 0, component_1, component_2) =
                    a(component_1, component_2);
            }
        }
        for (auto component_1 = 0; component_1 < 3; ++component_1) {
            for (auto component_2 = 0; component_2 < 3; ++component_2) {
                stiffness_matrix_terms(element, 0, 1, component_1, component_2) =
                    -a(component_1, component_2);
            }
        }
        for (auto component_1 = 0; component_1 < 3; ++component_1) {
            for (auto component_2 = 0; component_2 < 3; ++component_2) {
                stiffness_matrix_terms(element, 1, 0, component_1, component_2) =
                    -a(component_1, component_2);
            }
        }
        for (auto component_1 = 0; component_1 < 3; ++component_1) {
            for (auto component_2 = 0; component_2 < 3; ++component_2) {
                stiffness_matrix_terms(element, 1, 1, component_1, component_2) =
                    a(component_1, component_2);
            }
        }
    }
};

}  // namespace kynema_fmb::springs
