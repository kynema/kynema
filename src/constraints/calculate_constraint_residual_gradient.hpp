#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

#include "calculate_fixed_bc_3DOF_constraint.hpp"
#include "calculate_fixed_bc_constraint.hpp"
#include "calculate_prescribed_bc_3DOF_constraint.hpp"
#include "calculate_prescribed_bc_constraint.hpp"
#include "calculate_revolute_joint_constraint.hpp"
#include "calculate_revolute_joint_force.hpp"
#include "calculate_rigid_joint_3DOF_constraint.hpp"
#include "calculate_rigid_joint_constraint.hpp"
#include "calculate_rotation_control_constraint.hpp"
#include "constraint_type.hpp"

namespace kynema_fmb::constraints {

/**
 * @brief Top level kernel which calculates the residual and gradient contributions
 * of a constraint.
 *
 * @details This Kernel performs the action of identifying the constraint's type, loading
 * the required input variables into local variables, and then calling specialized
 * kernels that perform the actual calculations.
 */
template <typename DeviceType>
struct CalculateConstraintResidualGradient {
    using TransposeMatrix = KokkosBatched::SerialCopy<KokkosBatched::Trans::Transpose>;
    using Gemm = KokkosBatched::SerialGemm<
        KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Algo::Gemm::Default>;
    using CopyVector = KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose, 1>;
    using CopyMatrix = KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose, 2>;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<ConstraintType*> type_;
    ConstView<size_t*> base_node_index_;
    ConstView<size_t*> target_node_index_;
    ConstView<double* [3]> X0_;
    ConstView<double* [3][3]> axes_;
    ConstView<double* [7]> constraint_inputs_;
    ConstView<double* [6]> lambda_;
    ConstView<double* [6][6]> tangent_;
    ConstView<double* [7]> node_u_;
    View<double* [6]> res_;
    View<double* [6]> b_lambda_res_;
    View<double* [6]> t_lambda_res_;
    View<double* [6]> system_res_;
    View<double* [6][6]> b_grad_;
    View<double* [6][6]> t_grad_;
    View<double* [6][6]> b_grad_trans_;
    View<double* [6][6]> t_grad_trans_;

    KOKKOS_FUNCTION
    void FixedBC(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto t_grad_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculateFixedBCConstraint<DeviceType>::invoke(X0, t_node_u, res, t_grad);

        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void FixedBC3DOF(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto t_grad_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculateFixedBC3DOFConstraint<DeviceType>::invoke(X0, t_node_u, res, t_grad);

        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void PrescribedBC(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto inputs_data = Array<double, 7>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto t_grad_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto inputs = View<double[7]>(inputs_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(constraint_inputs_, constraint, ALL), inputs);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculatePrescribedBCConstraint<DeviceType>::invoke(X0, inputs, t_node_u, res, t_grad);

        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void PrescribedBC3DOF(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto inputs_data = Array<double, 7>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto t_grad_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto inputs = View<double[7]>(inputs_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(constraint_inputs_, constraint, ALL), inputs);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculatePrescribedBC3DOFConstraint<DeviceType>::invoke(X0, inputs, t_node_u, res, t_grad);

        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void RigidJoint(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto b_node_u_data = Array<double, 7>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto res_data = Array<double, 6>{};
        auto b_lambda_res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto b_grad_data = Array<double, 36>{};
        auto t_grad_data = Array<double, 36>{};
        auto b_grad_trans_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto base_tangent_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto b_grad_tan_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto b_node_u = View<double[7]>(b_node_u_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto b_lambda_res = View<double[6]>(b_lambda_res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto b_grad = View<double[6][6]>(b_grad_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto b_grad_trans = View<double[6][6]>(b_grad_trans_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto base_tangent = View<double[6][6]>(base_tangent_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto b_grad_tan = View<double[6][6]>(b_grad_tan_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto base_node_index = base_node_index_(constraint);
        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(node_u_, base_node_index, ALL), b_node_u);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(tangent_, base_node_index, ALL, ALL), base_tangent);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculateRigidJointConstraint<DeviceType>::invoke(
            X0, b_node_u, t_node_u, res, b_grad, t_grad
        );

        TransposeMatrix::invoke(b_grad, b_grad_trans);
        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., b_grad, base_tangent, 0., b_grad_tan);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., b_grad_trans, lambda, 0., b_lambda_res);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(b_lambda_res, subview(b_lambda_res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(b_grad_tan, subview(b_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(b_grad_trans, subview(b_grad_trans_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void RigidJoint3DOF(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto b_node_u_data = Array<double, 7>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto res_data = Array<double, 6>{};
        auto b_lambda_res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto b_grad_data = Array<double, 36>{};
        auto t_grad_data = Array<double, 36>{};
        auto b_grad_trans_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto base_tangent_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto b_grad_tan_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto b_node_u = View<double[7]>(b_node_u_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto b_lambda_res = View<double[6]>(b_lambda_res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto b_grad = View<double[6][6]>(b_grad_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto b_grad_trans = View<double[6][6]>(b_grad_trans_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto base_tangent = View<double[6][6]>(base_tangent_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto b_grad_tan = View<double[6][6]>(b_grad_tan_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto base_node_index = base_node_index_(constraint);
        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(node_u_, base_node_index, ALL), b_node_u);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(tangent_, base_node_index, ALL, ALL), base_tangent);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculateRigidJoint3DOFConstraint<DeviceType>::invoke(
            X0, b_node_u, t_node_u, res, b_grad, t_grad
        );

        TransposeMatrix::invoke(b_grad, b_grad_trans);
        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., b_grad, base_tangent, 0., b_grad_tan);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., b_grad_trans, lambda, 0., b_lambda_res);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(b_lambda_res, subview(b_lambda_res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(b_grad_tan, subview(b_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(b_grad_trans, subview(b_grad_trans_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void RevoluteJoint(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto inputs_data = Array<double, 7>{};
        auto b_node_u_data = Array<double, 7>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto axes_data = Array<double, 9>{};
        auto res_data = Array<double, 6>{};
        auto b_lambda_res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto system_res_data = Array<double, 6>{};
        auto b_grad_data = Array<double, 36>{};
        auto t_grad_data = Array<double, 36>{};
        auto b_grad_trans_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto base_tangent_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto b_grad_tan_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto inputs = View<double[7]>(inputs_data.data());
        const auto b_node_u = View<double[7]>(b_node_u_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto axes = View<double[3][3]>(axes_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto b_lambda_res = View<double[6]>(b_lambda_res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto system_res = View<double[6]>(system_res_data.data());
        const auto b_grad = View<double[6][6]>(b_grad_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto b_grad_trans = View<double[6][6]>(b_grad_trans_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto base_tangent = View<double[6][6]>(base_tangent_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto b_grad_tan = View<double[6][6]>(b_grad_tan_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto base_node_index = base_node_index_(constraint);
        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(constraint_inputs_, constraint, ALL), inputs);
        CopyVector::invoke(subview(node_u_, base_node_index, ALL), b_node_u);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(axes_, constraint, ALL, ALL), axes);
        CopyMatrix::invoke(subview(tangent_, base_node_index, ALL, ALL), base_tangent);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculateRevoluteJointConstraint<DeviceType>::invoke(
            X0, axes, b_node_u, t_node_u, res, b_grad, t_grad
        );
        CalculateRevoluteJointForce<DeviceType>::invoke(axes, inputs, t_node_u, system_res);

        TransposeMatrix::invoke(b_grad, b_grad_trans);
        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., b_grad, base_tangent, 0., b_grad_tan);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., b_grad_trans, lambda, 0., b_lambda_res);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(b_lambda_res, subview(b_lambda_res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyVector::invoke(system_res, subview(system_res_, constraint, ALL));
        CopyMatrix::invoke(b_grad_tan, subview(b_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(b_grad_trans, subview(b_grad_trans_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void RotationControl(size_t constraint) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using KokkosBlas::Experimental::serial_gemv;

        auto X0_data = Array<double, 3>{};
        auto inputs_data = Array<double, 7>{};
        auto b_node_u_data = Array<double, 7>{};
        auto t_node_u_data = Array<double, 7>{};
        auto lambda_data = Array<double, 6>{};
        auto axes_data = Array<double, 9>{};
        auto res_data = Array<double, 6>{};
        auto b_lambda_res_data = Array<double, 6>{};
        auto t_lambda_res_data = Array<double, 6>{};
        auto b_grad_data = Array<double, 36>{};
        auto t_grad_data = Array<double, 36>{};
        auto b_grad_trans_data = Array<double, 36>{};
        auto t_grad_trans_data = Array<double, 36>{};
        auto base_tangent_data = Array<double, 36>{};
        auto target_tangent_data = Array<double, 36>{};
        auto b_grad_tan_data = Array<double, 36>{};
        auto t_grad_tan_data = Array<double, 36>{};

        const auto X0 = View<double[3]>(X0_data.data());
        const auto inputs = View<double[7]>(inputs_data.data());
        const auto b_node_u = View<double[7]>(b_node_u_data.data());
        const auto t_node_u = View<double[7]>(t_node_u_data.data());
        const auto lambda = View<double[6]>(lambda_data.data());
        const auto axes = View<double[3][3]>(axes_data.data());
        const auto res = View<double[6]>(res_data.data());
        const auto b_lambda_res = View<double[6]>(b_lambda_res_data.data());
        const auto t_lambda_res = View<double[6]>(t_lambda_res_data.data());
        const auto b_grad = View<double[6][6]>(b_grad_data.data());
        const auto t_grad = View<double[6][6]>(t_grad_data.data());
        const auto b_grad_trans = View<double[6][6]>(b_grad_trans_data.data());
        const auto t_grad_trans = View<double[6][6]>(t_grad_trans_data.data());
        const auto base_tangent = View<double[6][6]>(base_tangent_data.data());
        const auto target_tangent = View<double[6][6]>(target_tangent_data.data());
        const auto b_grad_tan = View<double[6][6]>(b_grad_tan_data.data());
        const auto t_grad_tan = View<double[6][6]>(t_grad_tan_data.data());

        const auto base_node_index = base_node_index_(constraint);
        const auto target_node_index = target_node_index_(constraint);

        CopyVector::invoke(subview(X0_, constraint, ALL), X0);
        CopyVector::invoke(subview(constraint_inputs_, constraint, ALL), inputs);
        CopyVector::invoke(subview(node_u_, base_node_index, ALL), b_node_u);
        CopyVector::invoke(subview(node_u_, target_node_index, ALL), t_node_u);
        CopyVector::invoke(subview(lambda_, constraint, ALL), lambda);
        CopyMatrix::invoke(subview(axes_, constraint, ALL, ALL), axes);
        CopyMatrix::invoke(subview(tangent_, base_node_index, ALL, ALL), base_tangent);
        CopyMatrix::invoke(subview(tangent_, target_node_index, ALL, ALL), target_tangent);

        CalculateRotationControlConstraint<DeviceType>::invoke(
            X0, axes, inputs, b_node_u, t_node_u, res, b_grad, t_grad
        );

        TransposeMatrix::invoke(b_grad, b_grad_trans);
        TransposeMatrix::invoke(t_grad, t_grad_trans);
        Gemm::invoke(1., b_grad, base_tangent, 0., b_grad_tan);
        Gemm::invoke(1., t_grad, target_tangent, 0., t_grad_tan);
        serial_gemv('N', 1., b_grad_trans, lambda, 0., b_lambda_res);
        serial_gemv('N', 1., t_grad_trans, lambda, 0., t_lambda_res);

        CopyVector::invoke(res, subview(res_, constraint, ALL));
        CopyVector::invoke(b_lambda_res, subview(b_lambda_res_, constraint, ALL));
        CopyVector::invoke(t_lambda_res, subview(t_lambda_res_, constraint, ALL));
        CopyMatrix::invoke(b_grad_tan, subview(b_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_tan, subview(t_grad_, constraint, ALL, ALL));
        CopyMatrix::invoke(b_grad_trans, subview(b_grad_trans_, constraint, ALL, ALL));
        CopyMatrix::invoke(t_grad_trans, subview(t_grad_trans_, constraint, ALL, ALL));
    }

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        const auto constraint_type = type_(constraint);

        if (constraint_type == ConstraintType::FixedBC) {
            FixedBC(constraint);
        } else if (constraint_type == ConstraintType::FixedBC3DOFs) {
            FixedBC3DOF(constraint);
        } else if (constraint_type == ConstraintType::PrescribedBC) {
            PrescribedBC(constraint);
        } else if (constraint_type == ConstraintType::PrescribedBC3DOFs) {
            PrescribedBC3DOF(constraint);
        } else if (constraint_type == ConstraintType::RigidJoint) {
            RigidJoint(constraint);
        } else if (constraint_type == ConstraintType::RigidJoint6DOFsTo3DOFs) {
            RigidJoint3DOF(constraint);
        } else if (constraint_type == ConstraintType::RevoluteJoint) {
            RevoluteJoint(constraint);
        } else if (constraint_type == ConstraintType::RotationControl) {
            RotationControl(constraint);
        }
    };
};

}  // namespace kynema_fmb::constraints
