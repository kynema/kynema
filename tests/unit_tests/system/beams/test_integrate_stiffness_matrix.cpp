#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <gtest/gtest.h>

#include "system/beams/integrate_stiffness_matrix.hpp"
#include "test_calculate.hpp"

namespace {
void TestIntegrateStiffnessMatrix_1Element1Node1QP(
    const Kokkos::View<double[1][6][6]>::const_type& qp_Kuu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Puu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Cuu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Ouu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Quu,
    const std::array<double, 36>& exact_M_data, const double tolerance = 1e-12
) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_simd_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights = kynema::beams::tests::CreateView<double[number_of_qps]>(
        "qp_weights", std::array<double, 1>{2.}
    );
    const auto qp_jacobian = kynema::beams::tests::CreateView<double[number_of_qps]>(
        "qp_jacobian", std::array<double, 1>{3.}
    );
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::array<double, max_simd_size>{4.}
        );
    const auto shape_interp_deriv =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp_deriv", std::array<double, max_simd_size>{5.}
        );

    auto gbl_M = Kokkos::View<double[1][1][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateStiffnessMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_interp_deriv,
            .qp_Kuu_ = qp_Kuu,
            .qp_Puu_ = qp_Puu,
            .qp_Cuu_ = qp_Cuu,
            .qp_Ouu_ = qp_Ouu,
            .qp_Quu_ = qp_Quu,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    const auto exact_M =
        Kokkos::View<double[1][1][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M, tolerance);
}

constexpr std::array<double, 36> OneElement1Node1QP_Data() {
    return {0001., 0002., 0003., 0004., 0005., 0006., 1001., 1002., 1003., 1004., 1005., 1006.,
            2001., 2002., 2003., 2004., 2005., 2006., 3001., 3002., 3003., 3004., 3005., 3006.,
            4001., 4002., 4003., 4004., 4005., 4006., 5001., 5002., 5003., 5004., 5005., 5006.};
}

void TestIntegrateStiffnessMatrix_1Element1Node1QP_Kuu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Kuu", OneElement1Node1QP_Data()
    );
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data =
        std::array{000096., 000192., 000288., 000384., 000480., 000576., 096096., 096192., 096288.,
                   096384., 096480., 096576., 192096., 192192., 192288., 192384., 192480., 192576.,
                   288096., 288192., 288288., 288384., 288480., 288576., 384096., 384192., 384288.,
                   384384., 384480., 384576., 480096., 480192., 480288., 480384., 480480., 480576.};

    TestIntegrateStiffnessMatrix_1Element1Node1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element1Node1QP_Puu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Puu", OneElement1Node1QP_Data()
    );
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");

    constexpr auto exact_M_data =
        std::array{000040., 000080., 000120., 000160., 000200., 000240., 040040., 040080., 040120.,
                   040160., 040200., 040240., 080040., 080080., 080120., 080160., 080200., 080240.,
                   120040., 120080., 120120., 120160., 120200., 120240., 160040., 160080., 160120.,
                   160160., 160200., 160240., 200040., 200080., 200120., 200160., 200200., 200240.};

    TestIntegrateStiffnessMatrix_1Element1Node1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element1Node1QP_Quu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Quu", OneElement1Node1QP_Data()
    );
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data =
        std::array{000096., 000192., 000288., 000384., 000480., 000576., 096096., 096192., 096288.,
                   096384., 096480., 096576., 192096., 192192., 192288., 192384., 192480., 192576.,
                   288096., 288192., 288288., 288384., 288480., 288576., 384096., 384192., 384288.,
                   384384., 384480., 384576., 480096., 480192., 480288., 480384., 480480., 480576.};

    TestIntegrateStiffnessMatrix_1Element1Node1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

constexpr std::array<double, 36> OneElement1Node1QP_Cuu_Data() {
    return {03003., 03006., 03009., 03012., 03015., 03018., 06003., 06006., 06009.,
            06012., 06015., 06018., 09003., 09006., 09009., 09012., 09015., 09018.,
            12003., 12006., 12009., 12012., 12015., 12018., 15003., 15006., 15009.,
            15012., 15015., 15018., 18003., 18006., 18009., 18012., 18015., 18018.};
}

void TestIntegrateStiffnessMatrix_1Element1Node1QP_Cuu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Cuu", OneElement1Node1QP_Cuu_Data()
    );
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data =
        std::array{050050., 050100., 050150., 050200., 050250., 050300., 100050., 100100., 100150.,
                   100200., 100250., 100300., 150050., 150100., 150150., 150200., 150250., 150300.,
                   200050., 200100., 200150., 200200., 200250., 200300., 250050., 250100., 250150.,
                   250200., 250250., 250300., 300050., 300100., 300150., 300200., 300250., 300300.};

    TestIntegrateStiffnessMatrix_1Element1Node1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data, 1e-10
    );
}

constexpr std::array<double, 36> OneElement1Node1QP_Ouu_Data() {
    return {1001., 1002., 1003., 1004., 1005., 1006., 2001., 2002., 2003., 2004., 2005., 2006.,
            3001., 3002., 3003., 3004., 3005., 3006., 4001., 4002., 4003., 4004., 4005., 4006.,
            5001., 5002., 5003., 5004., 5005., 5006., 6001., 6002., 6003., 6004., 6005., 6006.};
}

void TestIntegrateStiffnessMatrix_1Element1Node1QP_Ouu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Ouu", OneElement1Node1QP_Ouu_Data()
    );

    constexpr auto exact_M_data =
        std::array{040040., 040080., 040120., 040160., 040200., 040240., 080040., 080080., 080120.,
                   080160., 080200., 080240., 120040., 120080., 120120., 120160., 120200., 120240.,
                   160040., 160080., 160120., 160160., 160200., 160240., 200040., 200080., 200120.,
                   200160., 200200., 200240., 240040., 240080., 240120., 240160., 240200., 240240.};

    TestIntegrateStiffnessMatrix_1Element1Node1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element2Nodes1QP(
    const Kokkos::View<double[1][6][6]>::const_type& qp_Kuu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Puu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Cuu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Ouu,
    const Kokkos::View<double[1][6][6]>::const_type& qp_Quu,
    const std::array<double, 144>& exact_M_data
) {
    constexpr auto number_of_nodes = size_t{2U};
    constexpr auto simd_width = Kokkos::Experimental::simd<double>::size();
    constexpr auto number_of_simd_nodes = (simd_width == 1) ? size_t{2U} : size_t{1U};
    constexpr auto number_of_qps = size_t{1U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights =
        kynema::beams::tests::CreateView<double[number_of_qps]>("qp_weights", std::array{1.});
    const auto qp_jacobian =
        kynema::beams::tests::CreateView<double[number_of_qps]>("qp_jacobian", std::array{1.});
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::vector<double>{1., 2., 0., 0., 0., 0., 0., 0.}
        );
    const auto shape_interp_deriv =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_deriv", std::vector<double>{1., 4., 0., 0., 0., 0., 0., 0.}
        );

    auto gbl_M = Kokkos::View<double[2][2][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateStiffnessMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_interp_deriv,
            .qp_Kuu_ = qp_Kuu,
            .qp_Puu_ = qp_Puu,
            .qp_Cuu_ = qp_Cuu,
            .qp_Ouu_ = qp_Ouu,
            .qp_Quu_ = qp_Quu,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    const auto exact_M =
        Kokkos::View<double[2][2][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

constexpr std::array<double, 36> OneElement2Nodes1QP_Data() {
    return {0001., 0002., 0003., 0004., 0005., 0006., 0101., 0102., 0103., 0104., 0105., 0106.,
            0201., 0202., 0203., 0204., 0205., 0206., 0301., 0302., 0303., 0304., 0305., 0306.,
            0401., 0402., 0403., 0404., 0405., 0406., 0501., 0502., 0503., 0504., 0505., 0506.};
}

void TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Puu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Puu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Puu", OneElement2Nodes1QP_Data()
    );
    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data = std::array<double, 144>{
        1.,    2.,    3.,    4.,    5.,    6.,    101.,  102.,  103.,  104.,  105.,  106.,
        201.,  202.,  203.,  204.,  205.,  206.,  301.,  302.,  303.,  304.,  305.,  306.,
        401.,  402.,  403.,  404.,  405.,  406.,  501.,  502.,  503.,  504.,  505.,  506.,
        4.,    8.,    12.,   16.,   20.,   24.,   404.,  408.,  412.,  416.,  420.,  424.,
        804.,  808.,  812.,  816.,  820.,  824.,  1204., 1208., 1212., 1216., 1220., 1224.,
        1604., 1608., 1612., 1616., 1620., 1624., 2004., 2008., 2012., 2016., 2020., 2024.,
        2.,    4.,    6.,    8.,    10.,   12.,   202.,  204.,  206.,  208.,  210.,  212.,
        402.,  404.,  406.,  408.,  410.,  412.,  602.,  604.,  606.,  608.,  610.,  612.,
        802.,  804.,  806.,  808.,  810.,  812.,  1002., 1004., 1006., 1008., 1010., 1012.,
        8.,    16.,   24.,   32.,   40.,   48.,   808.,  816.,  824.,  832.,  840.,  848.,
        1608., 1616., 1624., 1632., 1640., 1648., 2408., 2416., 2424., 2432., 2440., 2448.,
        3208., 3216., 3224., 3232., 3240., 3248., 4008., 4016., 4024., 4032., 4040., 4048.
    };

    TestIntegrateStiffnessMatrix_1Element2Nodes1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Quu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Quu", OneElement2Nodes1QP_Data()
    );
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

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

    TestIntegrateStiffnessMatrix_1Element2Nodes1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Cuu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Cuu", OneElement2Nodes1QP_Data()
    );
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data = std::array<double, 144>{
        1.,    2.,    3.,    4.,    5.,    6.,    101.,  102.,  103.,  104.,  105.,  106.,
        201.,  202.,  203.,  204.,  205.,  206.,  301.,  302.,  303.,  304.,  305.,  306.,
        401.,  402.,  403.,  404.,  405.,  406.,  501.,  502.,  503.,  504.,  505.,  506.,
        4.,    8.,    12.,   16.,   20.,   24.,   404.,  408.,  412.,  416.,  420.,  424.,
        804.,  808.,  812.,  816.,  820.,  824.,  1204., 1208., 1212., 1216., 1220., 1224.,
        1604., 1608., 1612., 1616., 1620., 1624., 2004., 2008., 2012., 2016., 2020., 2024.,
        4.,    8.,    12.,   16.,   20.,   24.,   404.,  408.,  412.,  416.,  420.,  424.,
        804.,  808.,  812.,  816.,  820.,  824.,  1204., 1208., 1212., 1216., 1220., 1224.,
        1604., 1608., 1612., 1616., 1620., 1624., 2004., 2008., 2012., 2016., 2020., 2024.,
        16.,   32.,   48.,   64.,   80.,   96.,   1616., 1632., 1648., 1664., 1680., 1696.,
        3216., 3232., 3248., 3264., 3280., 3296., 4816., 4832., 4848., 4864., 4880., 4896.,
        6416., 6432., 6448., 6464., 6480., 6496., 8016., 8032., 8048., 8064., 8080., 8096.
    };

    TestIntegrateStiffnessMatrix_1Element2Nodes1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Ouu() {
    constexpr auto number_of_qps = 1;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Puu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Ouu", OneElement2Nodes1QP_Data()
    );

    constexpr auto exact_M_data = std::array<double, 144>{
        1.,    2.,    3.,    4.,    5.,    6.,    101.,  102.,  103.,  104.,  105.,  106.,
        201.,  202.,  203.,  204.,  205.,  206.,  301.,  302.,  303.,  304.,  305.,  306.,
        401.,  402.,  403.,  404.,  405.,  406.,  501.,  502.,  503.,  504.,  505.,  506.,
        2.,    4.,    6.,    8.,    10.,   12.,   202.,  204.,  206.,  208.,  210.,  212.,
        402.,  404.,  406.,  408.,  410.,  412.,  602.,  604.,  606.,  608.,  610.,  612.,
        802.,  804.,  806.,  808.,  810.,  812.,  1002., 1004., 1006., 1008., 1010., 1012.,
        4.,    8.,    12.,   16.,   20.,   24.,   404.,  408.,  412.,  416.,  420.,  424.,
        804.,  808.,  812.,  816.,  820.,  824.,  1204., 1208., 1212., 1216., 1220., 1224.,
        1604., 1608., 1612., 1616., 1620., 1624., 2004., 2008., 2012., 2016., 2020., 2024.,
        8.,    16.,   24.,   32.,   40.,   48.,   808.,  816.,  824.,  832.,  840.,  848.,
        1608., 1616., 1624., 1632., 1640., 1648., 2408., 2416., 2424., 2432., 2440., 2448.,
        3208., 3216., 3224., 3232., 3240., 3248., 4008., 4016., 4024., 4032., 4040., 4048.
    };

    TestIntegrateStiffnessMatrix_1Element2Nodes1QP(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element1Node2QPs(
    const Kokkos::View<double[2][6][6]>::const_type& qp_Kuu,
    const Kokkos::View<double[2][6][6]>::const_type& qp_Puu,
    const Kokkos::View<double[2][6][6]>::const_type& qp_Cuu,
    const Kokkos::View<double[2][6][6]>::const_type& qp_Ouu,
    const Kokkos::View<double[2][6][6]>::const_type& qp_Quu,
    const std::array<double, 36>& exact_M_data
) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_simd_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{2U};
    constexpr auto max_simd_size = size_t{8U};

    const auto qp_weights =
        kynema::beams::tests::CreateView<double[number_of_qps]>("qp_weights", std::array{1., 3.});
    const auto qp_jacobian =
        kynema::beams::tests::CreateView<double[number_of_qps]>("qp_jacobian", std::array{2., 4.});
    const auto shape_interp =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp", std::array<double, max_simd_size * number_of_qps>{1., 1.}
        );
    const auto shape_interp_deriv =
        kynema::beams::tests::CreateLeftView<double[max_simd_size][number_of_qps]>(
            "shape_interp_deriv", std::array<double, max_simd_size * number_of_qps>{1., 1.}
        );

    auto gbl_M = Kokkos::View<double[1][1][6][6]>("global_M");

    const auto policy = Kokkos::RangePolicy(0, number_of_nodes * number_of_simd_nodes);
    const auto integrator =
        kynema::beams::IntegrateStiffnessMatrixElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = number_of_nodes,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_interp_deriv,
            .qp_Kuu_ = qp_Kuu,
            .qp_Puu_ = qp_Puu,
            .qp_Cuu_ = qp_Cuu,
            .qp_Ouu_ = qp_Ouu,
            .qp_Quu_ = qp_Quu,
            .gbl_M_ = gbl_M
        };
    Kokkos::parallel_for(policy, integrator);

    const auto exact_M =
        Kokkos::View<double[1][1][6][6], Kokkos::HostSpace>::const_type(exact_M_data.data());

    auto gbl_M_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gbl_M);
    kynema::beams::tests::CompareWithExpected(gbl_M_mirror, exact_M);
}

constexpr std::array<double, 72> OneElement1Node2QPs_Data() {
    return {10001., 10002., 10003., 10004., 10005., 10006., 10101., 10102., 10103., 10104., 10105.,
            10106., 10201., 10202., 10203., 10204., 10205., 10206., 10301., 10302., 10303., 10304.,
            10305., 10306., 10401., 10402., 10403., 10404., 10405., 10406., 10501., 10502., 10503.,
            10504., 10505., 10506., 20001., 20002., 20003., 20004., 20005., 20006., 20101., 20102.,
            20103., 20104., 20105., 20106., 20201., 20202., 20203., 20204., 20205., 20206., 20301.,
            20302., 20303., 20304., 20305., 20306., 20401., 20402., 20403., 20404., 20405., 20406.,
            20501., 20502., 20503., 20504., 20505., 20506.};
}

void TestIntegrateStiffnessMatrix_1Element1Node2QPs_Puu() {
    constexpr auto number_of_qps = 2;

    const auto qp_Puu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Puu", OneElement1Node2QPs_Data()
    );
    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Quu = Kokkos::View<double[number_of_qps][6][6]>("Quu");
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data =
        std::array{70004., 70008., 70012., 70016., 70020., 70024., 70404., 70408., 70412.,
                   70416., 70420., 70424., 70804., 70808., 70812., 70816., 70820., 70824.,
                   71204., 71208., 71212., 71216., 71220., 71224., 71604., 71608., 71612.,
                   71616., 71620., 71624., 72004., 72008., 72012., 72016., 72020., 72024.};

    TestIntegrateStiffnessMatrix_1Element1Node2QPs(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

void TestIntegrateStiffnessMatrix_1Element1Node2QPs_Quu() {
    constexpr auto number_of_qps = 2;

    const auto qp_Kuu = Kokkos::View<double[number_of_qps][6][6]>("Kuu");
    const auto qp_Puu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Quu = kynema::beams::tests::CreateView<double[number_of_qps][6][6]>(
        "qp_Quu", OneElement1Node2QPs_Data()
    );
    const auto qp_Cuu = Kokkos::View<double[number_of_qps][6][6]>("Cuu");
    const auto qp_Ouu = Kokkos::View<double[number_of_qps][6][6]>("Ouu");

    constexpr auto exact_M_data =
        std::array{260014., 260028., 260042., 260056., 260070., 260084., 261414., 261428., 261442.,
                   261456., 261470., 261484., 262814., 262828., 262842., 262856., 262870., 262884.,
                   264214., 264228., 264242., 264256., 264270., 264284., 265614., 265628., 265642.,
                   265656., 265670., 265684., 267014., 267028., 267042., 267056., 267070., 267084.};

    TestIntegrateStiffnessMatrix_1Element1Node2QPs(
        qp_Kuu, qp_Puu, qp_Cuu, qp_Ouu, qp_Quu, exact_M_data
    );
}

}  // namespace

namespace kynema::beams::tests {

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeOneQP_Kuu) {
    TestIntegrateStiffnessMatrix_1Element1Node1QP_Kuu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeOneQP_Puu) {
    TestIntegrateStiffnessMatrix_1Element1Node1QP_Puu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeOneQP_Quu) {
    TestIntegrateStiffnessMatrix_1Element1Node1QP_Quu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeOneQP_Cuu) {
    TestIntegrateStiffnessMatrix_1Element1Node1QP_Cuu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeOneQP_Ouu) {
    TestIntegrateStiffnessMatrix_1Element1Node1QP_Ouu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementTwoNodesOneQP_Puu) {
    TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Puu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementTwoNodesOneQP_Quu) {
    TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Quu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementTwoNodesOneQP_Cuu) {
    TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Cuu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementTwoNodesOneQP_Ouu) {
    TestIntegrateStiffnessMatrix_1Element2Nodes1QP_Ouu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeTwoQPs_Puu) {
    TestIntegrateStiffnessMatrix_1Element1Node2QPs_Puu();
}

TEST(IntegrateStiffnessMatrixTests, OneElementOneNodeTwoQPs_Quu) {
    TestIntegrateStiffnessMatrix_1Element1Node2QPs_Quu();
}

}  // namespace kynema::beams::tests
