#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Kernel for calculating the residual and system gradient for a rigid joint constraint
 * with three degrees of freedom
 */
template <typename DeviceType>
struct CalculateRigidJoint3DOFConstraint {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& X0, const ConstView<double[7]>& base_node_u,
        const ConstView<double[7]>& target_node_u, const View<double[6]>& residual_terms,
        const View<double[6][6]>& base_gradient_terms,
        const View<double[6][6]>& target_gradient_terms
    ) {
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;
        using CopyMatrix = KokkosBatched::SerialCopy<>;

        const auto u1_data = Array<double, 3>{base_node_u(0), base_node_u(1), base_node_u(2)};
        const auto R1_data =
            Array<double, 4>{base_node_u(3), base_node_u(4), base_node_u(5), base_node_u(6)};
        const auto u2_data = Array<double, 3>{target_node_u(0), target_node_u(1), target_node_u(2)};
        auto R1_X0_data = Array<double, 3>{};
        auto R1t_data = Array<double, 4>{};
        auto A_data = Array<double, 9>{};

        const auto u1 = ConstView<double[3]>{u1_data.data()};
        const auto R1 = ConstView<double[4]>{R1_data.data()};
        const auto u2 = ConstView<double[3]>{u2_data.data()};
        const auto R1_X0 = View<double[3]>{R1_X0_data.data()};
        const auto R1t = View<double[4]>{R1t_data.data()};
        const auto A = View<double[3][3]>{A_data.data()};

        //----------------------------------------------------------------------
        // Residual Vector
        //----------------------------------------------------------------------

        // Phi(0:3) = u2 + X0 - u1 - R1*X0
        math::QuaternionInverse(R1, R1t);
        math::RotateVectorByQuaternion(R1, X0, R1_X0);
        for (auto constraint = 0; constraint < 3; ++constraint) {
            residual_terms(constraint) =
                u2(constraint) + X0(constraint) - u1(constraint) - R1_X0(constraint);
        }

        //----------------------------------------------------------------------
        // Constraint Gradient Matrix
        //----------------------------------------------------------------------

        // Target Node gradients
        // B(0:3,0:3) = I
        for (auto constraint = 0; constraint < 3; ++constraint) {
            target_gradient_terms(constraint, constraint) = 1.;
        }

        // Base Node gradients
        // B(0:3,0:3) = -I
        for (auto constraint = 0; constraint < 3; ++constraint) {
            base_gradient_terms(constraint, constraint) = -1.;
        }

        // B(0:3,3:6) = tilde(R1*X0)
        math::VecTilde(R1_X0, A);
        CopyMatrix::invoke(A, subview(base_gradient_terms, make_pair(0, 3), make_pair(3, 6)));
    }
};
}  // namespace kynema_fmb::constraints
