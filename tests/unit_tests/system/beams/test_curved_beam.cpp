#include "test_curved_beam.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_SIMD.hpp>
#include <gtest/gtest.h>

#include "elements/beams/calculate_jacobian.hpp"
#include "math/interpolation.hpp"
#include "system/beams/calculate_Ouu.hpp"
#include "system/beams/calculate_Puu.hpp"
#include "system/beams/calculate_Quu.hpp"
#include "system/beams/calculate_force_FC.hpp"
#include "system/beams/calculate_force_FD.hpp"
#include "system/beams/integrate_residual_vector.hpp"
#include "system/beams/integrate_stiffness_matrix.hpp"
#include "system/masses/rotate_section_matrix.hpp"
#include "test_calculate.hpp"

/**
 * @file test_curved_beam.cpp
 * @brief Unit tests for S_t matrix, Residual, O, P, Q matrices, etc., of a beam element
 *
 * This file contains a series of unit tests that validate a curved beam problem consisting of
 * 3 nodes at GLL points (for second order polynomial). The validation dataset was created in
 * Mathematica during kynema theory development.
 *
 * The Mathematica notebooks used to generate the validation dataset are:
 * - Run1-RotationLibrary.nb
 * - Run2-ShapeFunctions.nb
 * - Run3-ReferenceLineDefinition.nb
 * - Run4-TestInterpolatedValues.nb
 * - Run5-DynamicIterationMatrix.nb
 *
 * The tests verify various aspects of system matrices construction for a Beam element such as:
 * - Lagrange polynomial/shape function -> LagrangePolynomialInterpWeights
 *   and derivative weights -> LagrangePolynomialDerivWeights
 * - Jacobian calculation -> CalculateJacobian
 * - Elastic forces calculation -> CalculateForceFc
 * - Damping forces calculation -> CalculateForceFd
 * - Residual vector integration -> beams::IntegrateResidualVectorElement
 * - Rotate sectional stiffness matrix -> beams::RotateSectionMatrix
 * - O matrix calculation -> CalculateOuu
 * - P matrix calculation -> CalculatePuu
 * - Q matrix calculation -> CalculateQuu
 * - Stiffness matrix integration -> beams::IntegrateStiffnessMatrixElement
 *
 * Ref: https://github.com/michaelasprague/KynemaTheory/tree/main/mathematica
 */

namespace {
void TestCalculateForceFc() {
    const auto Cuu =
        kynema::beams::tests::CreateView<double[6][6]>("Cuu", kynema::beams::tests::kCurvedBeamCuu);
    const auto strain = kynema::beams::tests::CreateView<double[6]>(
        "strain", kynema::beams::tests::kCurvedBeamStrain
    );

    const auto Fc = Kokkos::View<double[6]>("Fc");
    Kokkos::parallel_for(
        "CalculateForceFc", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateForceFC<Kokkos::DefaultExecutionSpace>::invoke(Cuu, strain, Fc);
        }
    );

    const auto Fc_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Fc);
    kynema::beams::tests::CompareWithExpected(Fc_mirror, kynema::beams::tests::kExpectedFc, 1e-6);
}

void TestCalculateForceFd() {
    const auto x0pupSS =
        kynema::beams::tests::CreateView<double[3][3]>("x0pupSS", kynema::beams::tests::kX0pupSS);
    const auto Fc =
        kynema::beams::tests::CreateView<double[6]>("Fc", kynema::beams::tests::kExpectedFc);

    const auto Fd = Kokkos::View<double[6]>("Fd");
    Kokkos::parallel_for(
        "CalculateForceFd", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateForceFD<Kokkos::DefaultExecutionSpace>::invoke(x0pupSS, Fc, Fd);
        }
    );

    const auto Fd_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Fd);
    kynema::beams::tests::CompareWithExpected(Fd_mirror, kynema::beams::tests::kExpectedFd, 1e-6);
}

void TestRotateSectionMatrixForCurvedBeam() {
    const auto xr =
        kynema::beams::tests::CreateView<double[4]>("xr", kynema::beams::tests::kCurvedBeamXr);
    const auto Cstar = kynema::beams::tests::CreateView<double[6][6]>(
        "Cstar", kynema::beams::tests::kCurvedBeamCstar
    );

    const auto Cuu = Kokkos::View<double[6][6]>("Cuu");
    Kokkos::parallel_for(
        "RotateSectionMatrix", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::masses::RotateSectionMatrix<Kokkos::DefaultExecutionSpace>::invoke(
                xr, Cstar, Cuu
            );
        }
    );

    const auto Cuu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Cuu);
    kynema::beams::tests::CompareWithExpected(Cuu_mirror, kynema::beams::tests::kExpectedCuu, 1e-6);
}

void TestCalculateOuu() {
    const auto Cuu =
        kynema::beams::tests::CreateView<double[6][6]>("Cuu", kynema::beams::tests::kCurvedBeamCuu);
    const auto x0pupSS =
        kynema::beams::tests::CreateView<double[3][3]>("x0pupSS", kynema::beams::tests::kX0pupSS);
    const auto M_tilde = kynema::beams::tests::CreateView<double[3][3]>(
        "M_tilde", kynema::beams::tests::kCurvedBeamM_tilde
    );
    const auto N_tilde = kynema::beams::tests::CreateView<double[3][3]>(
        "N_tilde", kynema::beams::tests::kCurvedBeamN_tilde
    );

    const auto Ouu = Kokkos::View<double[6][6]>("Ouu");
    Kokkos::parallel_for(
        "CalculateOuu", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateOuu<Kokkos::DefaultExecutionSpace>::invoke(
                Cuu, x0pupSS, M_tilde, N_tilde, Ouu
            );
        }
    );

    const auto Ouu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Ouu);
    kynema::beams::tests::CompareWithExpected(Ouu_mirror, kynema::beams::tests::kExpectedOuu, 1e-6);
}

void TestCalculatePuuForCurvedBeam() {
    const auto Cuu =
        kynema::beams::tests::CreateView<double[6][6]>("Cuu", kynema::beams::tests::kCurvedBeamCuu);
    const auto x0pupSS =
        kynema::beams::tests::CreateView<double[3][3]>("x0pupSS", kynema::beams::tests::kX0pupSS);
    const auto N_tilde = kynema::beams::tests::CreateView<double[3][3]>(
        "N_tilde", kynema::beams::tests::kCurvedBeamN_tilde
    );

    const auto Puu = Kokkos::View<double[6][6]>("Puu");
    Kokkos::parallel_for(
        "CalculatePuu", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculatePuu<Kokkos::DefaultExecutionSpace>::invoke(
                Cuu, x0pupSS, N_tilde, Puu
            );
        }
    );

    const auto Puu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Puu);
    kynema::beams::tests::CompareWithExpected(Puu_mirror, kynema::beams::tests::kExpectedPuu, 1e-6);
}

void TestCalculateQuuForCurvedBeam() {
    const auto Cuu =
        kynema::beams::tests::CreateView<double[6][6]>("Cuu", kynema::beams::tests::kCurvedBeamCuu);
    const auto x0pupSS =
        kynema::beams::tests::CreateView<double[3][3]>("x0pupSS", kynema::beams::tests::kX0pupSS);
    const auto N_tilde = kynema::beams::tests::CreateView<double[3][3]>(
        "N_tilde", kynema::beams::tests::kCurvedBeamN_tilde
    );

    const auto Quu = Kokkos::View<double[6][6]>("Quu");
    Kokkos::parallel_for(
        "CalculateQuu", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateQuu<Kokkos::DefaultExecutionSpace>::invoke(
                Cuu, x0pupSS, N_tilde, Quu
            );
        }
    );

    const auto Quu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Quu);
    kynema::beams::tests::CompareWithExpected(Quu_mirror, kynema::beams::tests::kExpectedQuu, 1e-6);
}

}  // namespace

namespace kynema::beams::tests {

TEST(CurvedBeamTests, LagrangePolynomialInterpWeight_SecondOrder_AtSpecifiedQPs) {
    std::vector<double> weights;
    // Test interpolation weights at each QP against expected data from Mathematica script
    for (auto qp : std::views::iota(0U, kNumQPs)) {
        math::LagrangePolynomialInterpWeights(kGaussQuadraturePoints[qp], kGLLNodes, weights);
        ASSERT_EQ(weights.size(), kNumNodes);
        EXPECT_NEAR(weights[0], kExpectedInterpWeights[qp][0], kDefaultTolerance);
        EXPECT_NEAR(weights[1], kExpectedInterpWeights[qp][1], kDefaultTolerance);
        EXPECT_NEAR(weights[2], kExpectedInterpWeights[qp][2], kDefaultTolerance);
    }
}

TEST(CurvedBeamTests, LagrangePolynomialDerivWeight_SecondOrder_AtSpecifiedQPs) {
    std::vector<double> deriv_weights;
    // Test derivative weights at each QP against expected data from Mathematica script
    for (auto qp : std::views::iota(0U, kNumQPs)) {
        math::LagrangePolynomialDerivWeights(kGaussQuadraturePoints[qp], kGLLNodes, deriv_weights);
        ASSERT_EQ(deriv_weights.size(), kNumNodes);
        EXPECT_NEAR(deriv_weights[0], kExpectedDerivWeights[qp][0], kDefaultTolerance);
        EXPECT_NEAR(deriv_weights[1], kExpectedDerivWeights[qp][1], kDefaultTolerance);
        EXPECT_NEAR(deriv_weights[2], kExpectedDerivWeights[qp][2], kDefaultTolerance);
    }
}

TEST(CurvedBeamTests, CalculateJacobianForCurvedBeam) {
    const auto node_position_rotation =
        CreateView<double[kNumElems][kNumNodes][7]>("node_position_rotation", kCurvedBeamNodes_data);
    const auto shape_derivative =
        CreateView<double[kNumElems][kNumNodes][kNumQPs]>("shape_derivative", kDerivWeightsFlat);
    const auto kNumNodes_per_elem =
        CreateView<size_t[kNumElems]>("number_of_nodes_per_element", std::array{kNumNodes});
    const auto kNumQPs_per_elem = CreateView<size_t[kNumElems]>(
        "number_of_quadrature_points_per_element", std::array{kNumQPs}
    );

    const auto qp_position_derivative =
        Kokkos::View<double[kNumElems][kNumQPs][3]>("position_derivative");
    const auto qp_jacobian = Kokkos::View<double[kNumElems][kNumQPs]>("jacobian");
    const auto calculate_jacobian = beams::CalculateJacobian<Kokkos::DefaultExecutionSpace>{
        .num_nodes_per_element = kNumNodes_per_elem,
        .num_qps_per_element = kNumQPs_per_elem,
        .shape_derivative = shape_derivative,
        .node_position_rotation = node_position_rotation,
        .qp_position_derivative = qp_position_derivative,
        .qp_jacobian = qp_jacobian
    };
    Kokkos::parallel_for("calculate_jacobian", 1, calculate_jacobian);
    const auto qp_jacobian_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_jacobian);

    // Validate the calculated jacobians against expected values from Mathematica script
    for (auto qp : std::views::iota(0U, kNumQPs)) {
        EXPECT_NEAR(qp_jacobian_mirror(0, qp), kExpectedJacobians[qp], kDefaultTolerance)
            << "Jacobian mismatch at quadrature point " << qp;
    }

    // Verify that position derivatives are unit vectors
    const auto qp_position_derivative_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_position_derivative);
    for (auto qp : std::views::iota(0U, kNumQPs)) {
        const auto magnitude = std::sqrt(
            (qp_position_derivative_mirror(0, qp, 0) * qp_position_derivative_mirror(0, qp, 0)) +
            (qp_position_derivative_mirror(0, qp, 1) * qp_position_derivative_mirror(0, qp, 1)) +
            (qp_position_derivative_mirror(0, qp, 2) * qp_position_derivative_mirror(0, qp, 2))
        );
        EXPECT_NEAR(magnitude, 1., kDefaultTolerance);  // unit vector
    }
}

TEST(CurvedBeamTests, CalculateForceFcForCurvedBeam) {
    TestCalculateForceFc();
}

TEST(CurvedBeamTests, CalculateForceFdForCurvedBeam) {
    TestCalculateForceFd();
}

TEST(CurvedBeamTests, IntegrateResidualVectorForCurvedBeam) {
    const auto qp_weights = CreateView<double[kNumQPs]>("qp_weights", kGaussQuadratureWeights);
    const auto qp_jacobian = CreateView<double[kNumQPs]>("qp_jacobian", kExpectedJacobians);
    const auto shape_interp =
        CreateLeftView<double[kNumNodes][kNumQPs]>("shape_interp", kInterpWeightsFlat);
    const auto shape_deriv =
        CreateLeftView<double[kNumNodes][kNumQPs]>("shape_deriv", kDerivWeightsFlat);

    const auto node_FX = Kokkos::View<double[kNumNodes][6]>("node_FX");
    const auto qp_Fc = CreateView<double[kNumQPs][6]>("qp_Fc", kFc);
    const auto qp_Fd = CreateView<double[kNumQPs][6]>("qp_Fd", kFd);
    const auto qp_Fi = CreateView<double[kNumQPs][6]>("qp_Fi", kFi);
    const auto qp_Fe = Kokkos::View<double[kNumQPs][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[kNumQPs][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[kNumElems][kNumNodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", kNumNodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = kNumQPs,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_deriv,
            .node_FX_ = node_FX,
            .qp_Fc_ = qp_Fc,
            .qp_Fd_ = qp_Fd,
            .qp_Fi_ = qp_Fi,
            .qp_Fe_ = qp_Fe,
            .qp_Fg_ = qp_Fg,
            .residual_vector_terms_ = residual_vector_terms
        }
    );

    const auto residual_vector_exact =
        Kokkos::View<double[kNumElems][kNumNodes][6], Kokkos::HostSpace>::const_type(
            kExpectedResidualVector.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, residual_vector_exact, 1e-8);
}

TEST(CurvedBeamTests, CalculateRotatedStiffnessMatrixForCurvedBeam) {
    TestRotateSectionMatrixForCurvedBeam();
}

TEST(CurvedBeamTests, CalculateOuuMatrixForCurvedBeam) {
    TestCalculateOuu();
}

TEST(CurvedBeamTests, CalculatePuuMatrixForCurvedBeam) {
    TestCalculatePuuForCurvedBeam();
}

TEST(CurvedBeamTests, CalculateQuuMatrixForCurvedBeam) {
    TestCalculateQuuForCurvedBeam();
}

TEST(CurvedBeamTests, IntegrateStiffnessMatrixForCurvedBeam) {
    constexpr auto max_simd_nodes = 8UL;
    const auto qp_weights = CreateView<double[kNumQPs]>("qp_weights", kGaussQuadratureWeights);
    const auto qp_jacobian = CreateView<double[kNumQPs]>("qp_jacobian", kExpectedJacobians);
    const auto shape_interp =
        CreateLeftView<double[max_simd_nodes][kNumQPs]>("shape_interp", kInterpWeightsFlat);
    const auto shape_deriv =
        CreateLeftView<double[max_simd_nodes][kNumQPs]>("shape_deriv", kDerivWeightsFlat);

    const auto qp_Kuu = CreateView<double[kNumQPs][6][6]>("qp_Kuu", kKuu);
    const auto qp_Puu = CreateView<double[kNumQPs][6][6]>("qp_Puu", kPuu);
    const auto qp_Cuu = CreateView<double[kNumQPs][6][6]>("qp_Cuu", kCuu);
    const auto qp_Ouu = CreateView<double[kNumQPs][6][6]>("qp_Ouu", kOuu);
    const auto qp_Quu = CreateView<double[kNumQPs][6][6]>("qp_Quu", kQuu);

    const auto stiffness_matrix_terms =
        Kokkos::View<double[kNumNodes][kNumNodes][6][6]>("stiffness_matrix_terms");

    constexpr auto simd_width = Kokkos::Experimental::simd<double>::size();
    constexpr auto extra_component = kNumNodes % simd_width == 0U ? 0U : 1U;
    constexpr auto simd_nodes = (kNumNodes / simd_width) + extra_component;
    const auto policy = Kokkos::RangePolicy(0, kNumNodes * simd_nodes);
    Kokkos::parallel_for(
        "IntegrateStiffnessMatrixElement", policy,
        beams::IntegrateStiffnessMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = kNumNodes,
            .num_qps = kNumQPs,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_deriv,
            .qp_Kuu_ = qp_Kuu,
            .qp_Puu_ = qp_Puu,
            .qp_Cuu_ = qp_Cuu,
            .qp_Ouu_ = qp_Ouu,
            .qp_Quu_ = qp_Quu,
            .gbl_M_ = stiffness_matrix_terms
        }
    );

    const auto stiffness_matrix_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), stiffness_matrix_terms);
    CompareWithExpected(stiffness_matrix_terms_mirror, kExpectedStiffnessMatrix, 1e-7);
}

}  // namespace kynema::beams::tests
