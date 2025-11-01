#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <gtest/gtest.h>

#include "system/beams/integrate_inertia_matrix.hpp"
#include "test_calculate.hpp"

namespace {

inline void IntegrateInertiaMatrix_TestOneElementOneNodeOneQP_Muu() {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_simd_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights = kynema::beams::tests::CreateView<double[number_of_qps]>(
        "weights", std::array<double, 1>{2.}
    );
    const auto qp_jacobian = kynema::beams::tests::CreateView<double[number_of_qps]>(
        "jacobian", std::array<double, 1>{3.}
    );
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::array<double, max_simd_size>{5.}
        );
    const auto qp_Muu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Muu", std::array{0001., 0002., 0003., 0004., 0005., 0006., 1001., 1002., 1003.,
                             1004., 1005., 1006., 2001., 2002., 2003., 2004., 2005., 2006.,
                             3001., 3002., 3003., 3004., 3005., 3006., 4001., 4002., 4003.,
                             4004., 4005., 4006., 5001., 5002., 5003., 5004., 5005., 5006.}
    );
    const auto qp_Guu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Guu", std::array{1001., 1002., 1003., 1004., 1005., 1006., 2001., 2002., 2003.,
                             2004., 2005., 2006., 4001., 4002., 4003., 4004., 4005., 4006.,
                             5001., 5002., 5003., 5004., 5005., 5006., 6001., 6002., 6003.,
                             6004., 6005., 6006., 7001., 7002., 7003., 7004., 7005., 7006.}
    );

    const auto gbl_M = Kokkos::View<double[1][1][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateInertiaMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .qp_Muu_ = qp_Muu,
            .qp_Guu_ = qp_Guu,
            .beta_prime_ = 1.,
            .gamma_prime_ = 0.,
            .gbl_M_ = gbl_M
        };

    Kokkos::parallel_for(policy, integrator);

    constexpr auto exact_M_data =
        std::array{000150., 000300., 000450., 000600., 000750., 000900., 150150., 150300., 150450.,
                   150600., 150750., 150900., 300150., 300300., 300450., 300600., 300750., 300900.,
                   450150., 450300., 450450., 450600., 450750., 450900., 600150., 600300., 600450.,
                   600600., 600750., 600900., 750150., 750300., 750450., 750600., 750750., 750900.};
    const auto exact_M =
        Kokkos::View<double[1][1][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

void IntegrateInertiaMatrix_TestOneElementOneNodeOneQP_Guu() {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_simd_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights =
        kynema::beams::tests::CreateView<double[number_of_qps]>("weights", std::array{2.});
    const auto qp_jacobian =
        kynema::beams::tests::CreateView<double[number_of_qps]>("jacobian", std::array{3.});
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::vector<double>{5., 0., 0., 0., 0., 0., 0., 0.}
        );
    const auto qp_Muu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Muu", std::array{0001., 0002., 0003., 0004., 0005., 0006., 1001., 1002., 1003.,
                             1004., 1005., 1006., 2001., 2002., 2003., 2004., 2005., 2006.,
                             3001., 3002., 3003., 3004., 3005., 3006., 4001., 4002., 4003.,
                             4004., 4005., 4006., 5001., 5002., 5003., 5004., 5005., 5006.}
    );
    const auto qp_Guu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Guu", std::array{1001., 1002., 1003., 1004., 1005., 1006., 2001., 2002., 2003.,
                             2004., 2005., 2006., 3001., 3002., 3003., 3004., 3005., 3006.,
                             4001., 4002., 4003., 4004., 4005., 4006., 5001., 5002., 5003.,
                             5004., 5005., 5006., 6001., 6002., 6003., 6004., 6005., 6006.}
    );

    const auto gbl_M = Kokkos::View<double[1][1][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateInertiaMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .qp_Muu_ = qp_Muu,
            .qp_Guu_ = qp_Guu,
            .beta_prime_ = 0.,
            .gamma_prime_ = 1.,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    constexpr auto exact_M_data =
        std::array{150150., 150300., 150450., 150600., 150750., 150900., 300150., 300300., 300450.,
                   300600., 300750., 300900., 450150., 450300., 450450., 450600., 450750., 450900.,
                   600150., 600300., 600450., 600600., 600750., 600900., 750150., 750300., 750450.,
                   750600., 750750., 750900., 900150., 900300., 900450., 900600., 900750., 900900.};
    const auto exact_M =
        Kokkos::View<double[1][1][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    Kokkos::deep_copy(gbl_M_mirror, gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

void IntegrateInertiaMatrix_TestOneElementTwoNodesOneQP() {
    constexpr auto number_of_nodes = size_t{2U};
    constexpr auto simd_width = Kokkos::Experimental::simd<double>::size();
    constexpr auto number_of_simd_nodes = (simd_width == 1) ? size_t{2U} : size_t{1U};
    constexpr auto number_of_qps = size_t{1U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights =
        kynema::beams::tests::CreateView<double[number_of_qps]>("weights", std::array{1.});
    const auto qp_jacobian =
        kynema::beams::tests::CreateView<double[number_of_qps]>("jacobian", std::array{1.});
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp",
            std::vector<double>{1., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}
        );
    const auto qp_Muu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Muu", std::array{001., 002., 003., 004., 005., 006., 101., 102., 103., 104., 105., 106.,
                             201., 202., 203., 204., 205., 206., 301., 302., 303., 304., 305., 306.,
                             401., 402., 403., 404., 405., 406., 501., 502., 503., 504., 505., 506.}
    );
    const auto qp_Guu = Kokkos::View<double[number_of_qps][6][6]>("Guu");

    const auto gbl_M = Kokkos::View<double[2][2][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateInertiaMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .qp_Muu_ = qp_Muu,
            .qp_Guu_ = qp_Guu,
            .beta_prime_ = 1.,
            .gamma_prime_ = 0.,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    constexpr auto exact_M_data = std::array<double, 144>{
        1.,    2.,    3.,    4.,    5.,    6.,    101.,  102.,  103.,  104.,  105.,  106.,
        201.,  202.,  203.,  204.,  205.,  206.,  301.,  302.,  303.,  304.,  305.,  306.,
        401.,  402.,  403.,  404.,  405.,  406.,  501.,  502.,  503.,  504.,  505.,  506.,
        2.,    4.,    6.,    8.,    10.,   12.,   202.,  204.,  206.,  208.,  210.,  212.,
        402.,  404.,  406.,  408.,  410.,  412.,  602.,  604.,  606.,  608.,  610.,  612.,
        802.,  804.,  806.,  808.,  810.,  812.,  1002., 1004., 1006., 1008., 1010., 1012.,
        2.,    4.,    6.,    8.,    10.,   12.,   202.,  204.,  206.,  208.,  210.,  212.,
        402.,  404.,  406.,  408.,  410.,  412.,  602.,  604.,  606.,  608.,  610.,  612.,
        802.,  804.,  806.,  808.,  810.,  812.,  1002., 1004., 1006., 1008., 1010., 1012.,
        4.,    8.,    12.,   16.,   20.,   24.,   404.,  408.,  412.,  416.,  420.,  424.,
        804.,  808.,  812.,  816.,  820.,  824.,  1204., 1208., 1212., 1216., 1220., 1224.,
        1604., 1608., 1612., 1616., 1620., 1624., 2004., 2008., 2012., 2016., 2020., 2024.
    };

    const auto exact_M =
        Kokkos::View<double[2][2][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

void IntegrateInertiaMatrix_TestOneElementOneNodeTwoQPs() {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_simd_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{2U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights =
        kynema::beams::tests::CreateView<double[number_of_qps]>("weights", std::array{9., 1.});
    const auto qp_jacobian =
        kynema::beams::tests::CreateView<double[number_of_qps]>("jacobian", std::array{1., 4.});
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::array<double, max_simd_size * number_of_qps>{1. / 3., 1. / 2.}
        );
    const auto qp_Muu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Muu", std::array{00001., 00002., 00003., 00004., 00005., 00006., 00011., 00012., 00013.,
                             00014., 00015., 00016., 00021., 00022., 00023., 00024., 00025., 00026.,
                             00031., 00032., 00033., 00034., 00035., 00036., 00041., 00042., 00043.,
                             00044., 00045., 00046., 00051., 00052., 00053., 00054., 00055., 00056.,
                             01000., 02000., 03000., 04000., 05000., 06000., 11000., 12000., 13000.,
                             14000., 15000., 16000., 21000., 22000., 23000., 24000., 25000., 26000.,
                             31000., 32000., 33000., 34000., 35000., 36000., 41000., 42000., 43000.,
                             44000., 45000., 46000., 51000., 52000., 53000., 54000., 55000., 56000.}
    );
    const auto qp_Guu = Kokkos::View<double[number_of_qps][6][6]>("Guu");

    const auto gbl_M = Kokkos::View<double[1][1][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateInertiaMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .qp_Muu_ = qp_Muu,
            .qp_Guu_ = qp_Guu,
            .beta_prime_ = 1.,
            .gamma_prime_ = 0.,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    constexpr auto exact_M_data =
        std::array{01001., 02002., 03003., 04004., 05005., 06006., 11011., 12012., 13013.,
                   14014., 15015., 16016., 21021., 22022., 23023., 24024., 25025., 26026.,
                   31031., 32032., 33033., 34034., 35035., 36036., 41041., 42042., 43043.,
                   44044., 45045., 46046., 51051., 52052., 53053., 54054., 55055., 56056.};
    const auto exact_M =
        Kokkos::View<double[1][1][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

void IntegrateInertiaMatrix_TestOneElementOneNodeOneQP_WithMultiplicationFactor() {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_simd_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights =
        kynema::beams::tests::CreateView<double[number_of_qps]>("weights", std::array{1.});
    const auto qp_jacobian =
        kynema::beams::tests::CreateView<double[number_of_qps]>("jacobian", std::array{1.});
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::vector<double>{1., 0., 0., 0., 0., 0., 0., 0.}
        );
    const auto qp_Muu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Muu", std::array{0001., 0002., 0003., 0004., 0005., 0006., 1001., 1002., 1003.,
                             1004., 1005., 1006., 2001., 2002., 2003., 2004., 2005., 2006.,
                             3001., 3002., 3003., 3004., 3005., 3006., 4001., 4002., 4003.,
                             4004., 4005., 4006., 5001., 5002., 5003., 5004., 5005., 5006.}
    );
    const auto qp_Guu = Kokkos::View<double[number_of_qps][6][6]>("Guu");

    const auto multiplication_factor = 5.;

    const auto gbl_M = Kokkos::View<double[1][1][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateInertiaMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .qp_Muu_ = qp_Muu,
            .qp_Guu_ = qp_Guu,
            .beta_prime_ = multiplication_factor,
            .gamma_prime_ = 0.,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    constexpr auto exact_M_data =
        std::array{00005., 00010., 00015., 00020., 00025., 00030., 05005., 05010., 05015.,
                   05020., 05025., 05030., 10005., 10010., 10015., 10020., 10025., 10030.,
                   15005., 15010., 15015., 15020., 15025., 15030., 20005., 20010., 20015.,
                   20020., 20025., 20030., 25005., 25010., 25015., 25020., 25025., 25030.};
    const auto exact_M =
        Kokkos::View<double[1][1][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

}  // namespace
namespace kynema::beams::tests {

TEST(IntegrateInertiaMatrixTests, OneElementOneNodeOneQP_Muu) {
    IntegrateInertiaMatrix_TestOneElementOneNodeOneQP_Muu();
}

TEST(IntegrateInertiaMatrixTests, OneElementOneNodeOneQP_Guu) {
    IntegrateInertiaMatrix_TestOneElementOneNodeOneQP_Guu();
}

TEST(IntegrateInertiaMatrixTests, OneElementTwoNodesOneQP) {
    IntegrateInertiaMatrix_TestOneElementTwoNodesOneQP();
}

TEST(IntegrateInertiaMatrixTests, OneElementOneNodeTwoQPs) {
    IntegrateInertiaMatrix_TestOneElementOneNodeTwoQPs();
}

TEST(IntegrateInertiaMatrixTests, OneElementOneNodeOneQP_WithMultiplicationFactor_TeamPolicy) {
    IntegrateInertiaMatrix_TestOneElementOneNodeOneQP_WithMultiplicationFactor();
}

}  // namespace kynema::beams::tests
