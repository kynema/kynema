#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "math/matrix_operations.hpp"
#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief A Kernel which calculates the residual and gradient contributions of a rotation
 * control constraint
 */
template <typename DeviceType>
struct CalculateRotationControlConstraint {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& X0, const ConstView<double[3][3]>& axes,
        const ConstView<double[7]>& inputs, const ConstView<double[7]>& base_node_u,
        const ConstView<double[7]>& target_node_u, const View<double[6]>& residual_terms,
        const View<double[6][6]>& base_gradient_terms,
        const View<double[6][6]>& target_gradient_terms
    ) {
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;
        using CopyVector = KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose, 1>;
        using CopyMatrixTranspose = KokkosBatched::SerialCopy<KokkosBatched::Trans::Transpose>;
        using CopyMatrix = KokkosBatched::SerialCopy<>;

        const auto ub_data = Array<double, 3>{base_node_u(0), base_node_u(1), base_node_u(2)};
        const auto Rb_data =
            Array<double, 4>{base_node_u(3), base_node_u(4), base_node_u(5), base_node_u(6)};
        const auto ut_data = Array<double, 3>{target_node_u(0), target_node_u(1), target_node_u(2)};
        const auto Rt_data =
            Array<double, 4>{target_node_u(3), target_node_u(4), target_node_u(5), target_node_u(6)};
        auto AX_data = Array<double, 3>{};
        auto RV_data = Array<double, 3>{};
        auto Rc_data = Array<double, 4>{};
        auto RcT_data = Array<double, 4>{};
        auto RbT_data = Array<double, 4>{};
        auto Rb_X0_data = Array<double, 4>{};
        auto Rt_RcT_data = Array<double, 4>{};
        auto Rt_RcT_RbT_data = Array<double, 4>{};
        auto A_data = Array<double, 9>{};
        auto C_data = Array<double, 9>{};
        auto V3_data = Array<double, 3>{};

        const auto ub = ConstView<double[3]>{ub_data.data()};
        const auto Rb = ConstView<double[4]>{Rb_data.data()};
        const auto ut = ConstView<double[3]>{ut_data.data()};
        const auto Rt = ConstView<double[4]>{Rt_data.data()};
        const auto AX = View<double[3]>{AX_data.data()};
        const auto RV = View<double[3]>{RV_data.data()};
        const auto Rc = View<double[4]>{Rc_data.data()};
        const auto RcT = View<double[4]>{RcT_data.data()};
        const auto RbT = View<double[4]>{RbT_data.data()};
        const auto Rb_X0 = View<double[4]>{Rb_X0_data.data()};
        const auto Rt_RcT = View<double[4]>{Rt_RcT_data.data()};
        const auto Rt_RcT_RbT = View<double[4]>{Rt_RcT_RbT_data.data()};
        const auto A = View<double[3][3]>{A_data.data()};
        const auto C = View<double[3][3]>{C_data.data()};
        const auto V3 = View<double[3]>{V3_data.data()};

        //----------------------------------------------------------------------
        // Position residual
        //----------------------------------------------------------------------

        // Phi(0:3) = ut + X0 - ub - Rb*X0
        math::RotateVectorByQuaternion(Rb, X0, Rb_X0);
        for (auto i = 0; i < 3; ++i) {
            residual_terms(i) = ut(i) + X0(i) - ub(i) - Rb_X0(i);
        }

        //----------------------------------------------------------------------
        // Rotation residual
        //----------------------------------------------------------------------

        auto rotation_command = inputs(0);

        // Copy rotation axis for this constraint
        for (auto i = 0U; i < 3U; ++i) {
            AX(i) = axes(0, i);
            RV(i) = AX(i) * rotation_command;
        }

        // Convert scaled axis to quaternion and calculate inverse
        math::RotationVectorToQuaternion(RV, Rc);
        math::QuaternionInverse(Rc, RcT);

        // Phi(3:6) = axial(Rt*inv(Rc)*inv(Rb))
        math::QuaternionInverse(Rb, RbT);
        math::QuaternionCompose(Rt, RcT, Rt_RcT);
        math::QuaternionCompose(Rt_RcT, RbT, Rt_RcT_RbT);
        math::QuaternionToRotationMatrix(Rt_RcT_RbT, C);
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

        // B(3:6,3:6) = AX(Rb*Rc*inv(Rt)) = transpose(AX(Rt*inv(Rc)*inv(Rb)))
        math::AX_Matrix(C, A);
        CopyMatrixTranspose::invoke(
            A, subview(target_gradient_terms, make_pair(3, 6), make_pair(3, 6))
        );

        //---------------------------------
        // Base Node
        //---------------------------------
        // B(0:3,0:3) = -I
        for (auto component = 0; component < 3; ++component) {
            base_gradient_terms(component, component) = -1.;
        }

        // B(0:3,3:6) = tilde(Rb*X0)
        math::VecTilde(Rb_X0, A);
        CopyMatrix::invoke(A, subview(base_gradient_terms, make_pair(0, 3), make_pair(3, 6)));

        // B(3:6,3:6) = -AX(Rt*inv(Rc)*inv(Rb))
        math::AX_Matrix(C, A);
        for (auto component_1 = 0; component_1 < 3; ++component_1) {
            for (auto component_2 = 0; component_2 < 3; ++component_2) {
                base_gradient_terms(component_1 + 3, component_2 + 3) = -A(component_1, component_2);
            }
        }
    }
};
}  // namespace kynema_fmb::constraints
