#pragma once

#include <Kokkos_Core.hpp>

#include "calculate_dissipation_quadrature_point_values.hpp"
#include "calculate_inertial_quadrature_point_values.hpp"
#include "calculate_stiffness_quadrature_point_values.hpp"
#include "calculate_system_matrix.hpp"
#include "integrate_inertia_matrix.hpp"
#include "integrate_residual_vector.hpp"
#include "integrate_stiffness_matrix.hpp"
#include "update_node_state.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateQuadraturePointValues {
    using TeamPolicy = typename Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    double beta_prime_;
    double gamma_prime_;
    ConstView<double* [7]> Q;
    ConstView<double* [6]> V;
    ConstView<double* [6]> A;
    ConstView<double* [6][6]> tangent;
    ConstView<size_t**> node_state_indices;
    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t*> num_qps_per_element;
    ConstView<double**> qp_weight_;
    ConstView<double**> qp_jacobian_;
    ConstView<double***> shape_interp_;
    ConstView<double***> shape_deriv_;
    ConstView<double[3]> gravity_;
    ConstView<double** [6]> node_FX_;
    ConstView<double** [4]> qp_r0_;
    ConstView<double** [3]> qp_x0_;
    ConstView<double** [3]> qp_x0_prime_;
    ConstView<double** [6][6]> qp_Mstar_;
    ConstView<double** [6][6]> qp_Cstar_;
    View<double** [6]> qp_FE_;
    View<double** [6]> residual_vector_terms_;
    View<double*** [6][6]> system_matrix_terms_;

    KOKKOS_FUNCTION
    void operator()(member_type member) const {
        using simd_type = Kokkos::Experimental::simd<double>;
        using Kokkos::ALL;
        using Kokkos::make_pair;
        using Kokkos::parallel_for;
        using Kokkos::subview;
        using Kokkos::TeamVectorRange;
        using CopyMatrix = KokkosBatched::TeamVectorCopy<member_type>;
        using CopyVector =
            KokkosBatched::TeamVectorCopy<member_type, KokkosBatched::Trans::NoTranspose, 1>;

        const auto element = static_cast<size_t>(member.league_rank());
        const auto num_nodes = num_nodes_per_element(element);
        const auto num_qps = num_qps_per_element(element);
        constexpr auto width = simd_type::size();
        const auto extra_component = num_nodes % width == 0U ? 0U : 1U;
        const auto simd_nodes = (num_nodes / width) + extra_component;
        const auto padded_num_nodes = simd_nodes * width;

        const auto qp_pair = make_pair(0UL, num_qps);
        const auto node_pair = make_pair(0UL, num_nodes);
        const auto qp_range = TeamVectorRange(member, num_qps);
        const auto node_range = TeamVectorRange(member, num_nodes);
        const auto node_squared_range = TeamVectorRange(member, num_nodes * num_nodes);
        const auto node_squared_simd_range = TeamVectorRange(member, num_nodes * simd_nodes);

        const auto shape_interp =
            LeftView<double**>(member.team_scratch(0), padded_num_nodes, num_qps);
        const auto shape_deriv =
            LeftView<double**>(member.team_scratch(0), padded_num_nodes, num_qps);

        const auto qp_weight = View<double*>(member.team_scratch(0), num_qps);
        const auto qp_jacobian = View<double*>(member.team_scratch(0), num_qps);

        const auto node_u = View<double* [7]>(member.team_scratch(1), num_nodes);
        const auto node_u_dot = View<double* [6]>(member.team_scratch(1), num_nodes);
        const auto node_u_ddot = View<double* [6]>(member.team_scratch(1), num_nodes);
        const auto node_FX = View<double* [6]>(member.team_scratch(1), num_nodes);
        const auto qp_Fc = View<double* [6]>(member.team_scratch(1), num_qps);
        const auto qp_Fd = View<double* [6]>(member.team_scratch(1), num_qps);
        const auto qp_Fi = View<double* [6]>(member.team_scratch(1), num_qps);
        const auto qp_Fe = View<double* [6]>(member.team_scratch(1), num_qps);
        const auto qp_Fg = View<double* [6]>(member.team_scratch(1), num_qps);
        const auto qp_FD1 = View<double* [6]>(member.team_scratch(1), num_qps);
        const auto qp_FD2 = View<double* [6]>(member.team_scratch(1), num_qps);

        const auto qp_Kuu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_Puu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_Cuu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_Ouu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_Quu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_Muu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_Guu = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_DD1 = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_DD2 = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_GD1 = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_GD2 = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_PD2 = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_KD1 = View<double* [6][6]>(member.team_scratch(1), num_qps);
        const auto qp_KD2 = View<double* [6][6]>(member.team_scratch(1), num_qps);

        const auto stiffness_matrix_terms =
            View<double** [6][6]>(member.team_scratch(1), num_nodes, num_nodes);
        const auto inertia_matrix_terms =
            View<double** [6][6]>(member.team_scratch(1), num_nodes, num_nodes);

        CopyMatrix::invoke(
            member, subview(shape_interp_, element, node_pair, qp_pair),
            subview(shape_interp, node_pair, qp_pair)
        );
        CopyMatrix::invoke(
            member, subview(shape_deriv_, element, node_pair, qp_pair),
            subview(shape_deriv, node_pair, qp_pair)
        );
        CopyMatrix::invoke(member, subview(qp_FE_, element, qp_pair, ALL), qp_Fe);
        CopyMatrix::invoke(member, subview(node_FX_, element, node_pair, ALL), node_FX);

        CopyVector::invoke(member, subview(qp_weight_, element, qp_pair), qp_weight);
        CopyVector::invoke(member, subview(qp_jacobian_, element, qp_pair), qp_jacobian);

        const auto node_state_updater = beams::UpdateNodeStateElement<DeviceType>{
            element, node_state_indices, node_u, node_u_dot, node_u_ddot, Q, V, A
        };
        parallel_for(node_range, node_state_updater);
        member.team_barrier();

        const auto inertia_quad_point_calculator =
            beams::CalculateInertialQuadraturePointValues<DeviceType>{
                element,     shape_interp, gravity_, qp_r0_, qp_Mstar_, node_u, node_u_dot,
                node_u_ddot, qp_Fi,        qp_Fg,    qp_Muu, qp_Guu,    qp_Kuu
            };
        parallel_for(qp_range, inertia_quad_point_calculator);

        const auto stiffness_quad_point_calculator =
            beams::CalculateStiffnessQuadraturePointValues<DeviceType>{
                element, qp_jacobian, shape_interp, shape_deriv, qp_r0_, qp_x0_prime_, qp_Cstar_,
                node_u,  qp_Fc,       qp_Fd,        qp_Cuu,      qp_Ouu, qp_Puu,       qp_Quu
            };
        parallel_for(qp_range, stiffness_quad_point_calculator);

        /*
            const auto dissipation_quad_point_calculator =
                beams::CalculateDissipationQuadraturePointValues<DeviceType>{
                    element, qp_jacobian, shape_interp, shape_deriv, qp_r0_, qp_FD1, qp_FD2, qp_DD1,
           qp_DD2, qp_GD1, qp_GD2, qp_PD2, qp_KD1, qp_KD2
                };
            parallel_for(qp_range, dissipation_quad_point_calculator);
          */
        member.team_barrier();

        const auto residual_integrator = beams::IntegrateResidualVectorElement<DeviceType>{
            element, num_qps, qp_weight, qp_jacobian, shape_interp, shape_deriv,           node_FX,
            qp_Fc,   qp_Fd,   qp_Fi,     qp_Fe,       qp_Fg,        residual_vector_terms_
        };
        parallel_for(node_range, residual_integrator);

        const auto stiffness_matrix_integrator = beams::IntegrateStiffnessMatrixElement<DeviceType>{
            element, num_nodes, num_qps, qp_weight, qp_jacobian, shape_interp,          shape_deriv,
            qp_Kuu,  qp_Puu,    qp_Cuu,  qp_Ouu,    qp_Quu,      stiffness_matrix_terms
        };
        parallel_for(node_squared_simd_range, stiffness_matrix_integrator);

        const auto inertia_matrix_integrator = beams::IntegrateInertiaMatrixElement<DeviceType>{
            element, num_nodes, num_qps,     qp_weight,    qp_jacobian,         shape_interp,
            qp_Muu,  qp_Guu,    beta_prime_, gamma_prime_, inertia_matrix_terms
        };
        parallel_for(node_squared_simd_range, inertia_matrix_integrator);
        member.team_barrier();

        const auto system_matrix_calculator = beams::CalculateSystemMatrix<DeviceType>{
            element,
            num_nodes,
            tangent,
            node_state_indices,
            stiffness_matrix_terms,
            inertia_matrix_terms,
            system_matrix_terms_
        };
        parallel_for(node_squared_range, system_matrix_calculator);
    }
};

}  // namespace kynema::beams
