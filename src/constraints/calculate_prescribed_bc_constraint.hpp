#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "math/matrix_operations.hpp"
#include "math/quaternion_operations.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Kernel for calculating the residual and system gradient for a Prescribed BC constraint
 * with six degrees of freedom
 */
template <typename DeviceType>
struct CalculatePrescribedBCConstraint {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& X0, const ConstView<double[7]>& inputs,
        const ConstView<double[7]>& node_u, const View<double[6]>& residual_terms,
        const View<double[6][6]>& target_gradient_terms
    ) {
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;
        using CopyVector = KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose, 1>;
        using CopyTransposeMatrix = KokkosBatched::SerialCopy<KokkosBatched::Trans::Transpose>;

        const auto u1_data = Array<double, 3>{inputs(0), inputs(1), inputs(2)};
        const auto R1_data = Array<double, 4>{inputs(3), inputs(4), inputs(5), inputs(6)};
        const auto u2_data = Array<double, 3>{node_u(0), node_u(1), node_u(2)};
        const auto R2_data = Array<double, 4>{node_u(3), node_u(4), node_u(5), node_u(6)};
        auto R1t_data = Array<double, 4>{};
        auto R1_X0_data = Array<double, 4>{};
        auto R2_R1t_data = Array<double, 4>{};
        auto A_data = Array<double, 9>{};
        auto C_data = Array<double, 9>{};
        auto V3_data = Array<double, 3>{};

        const auto u1 = ConstView<double[3]>{u1_data.data()};
        const auto R1 = ConstView<double[4]>{R1_data.data()};
        const auto u2 = ConstView<double[3]>{u2_data.data()};
        const auto R2 = ConstView<double[4]>{R2_data.data()};
        const auto R1t = View<double[4]>{R1t_data.data()};
        const auto R1_X0 = View<double[4]>{R1_X0_data.data()};
        const auto R2_R1t = View<double[4]>{R2_R1t_data.data()};
        const auto A = View<double[3][3]>{A_data.data()};
        const auto C = View<double[3][3]>{C_data.data()};
        const auto V3 = View<double[3]>{V3_data.data()};

        //----------------------------------------------------------------------
        // Residual Vector
        //----------------------------------------------------------------------

        // Phi(0:3) = u2 + X0 - u1 - R1*X0
        math::QuaternionInverse(R1, R1t);
        math::RotateVectorByQuaternion(R1, X0, R1_X0);
        for (auto component = 0; component < 3; ++component) {
            residual_terms(component) =
                u2(component) + X0(component) - u1(component) - R1_X0(component);
        }

        // Angular residual
        // Phi(3:6) = axial(R2*inv(RC)*inv(R1))
        math::QuaternionCompose(R2, R1t, R2_R1t);
        math::QuaternionToRotationMatrix(R2_R1t, C);
        math::AxialVectorOfMatrix(C, V3);
        CopyVector::invoke(V3, subview(residual_terms, make_pair(3, 6)));

        //----------------------------------------------------------------------
        // Constraint Gradient Matrix
        //----------------------------------------------------------------------

        //---------------------------------
        // Target Node
        //---------------------------------

        // B(0:3,0:3) = I
        for (auto component = 0; component < 3; ++component) {
            target_gradient_terms(component, component) = 1.;
        }

        // B(3:6,3:6) = AX(R1*RC*inv(R2)) = transpose(AX(R2*inv(RC)*inv(R1)))
        math::AX_Matrix(C, A);
        CopyTransposeMatrix::invoke(
            A, subview(target_gradient_terms, make_pair(3, 6), make_pair(3, 6))
        );
    }
};
}  // namespace kynema_fmb::constraints
