#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "system/beams/integrate_residual_vector.hpp"
#include "test_calculate.hpp"

namespace kynema::beams::tests {

TEST(IntegrateResidualVector, OneElementOneNodeOneQP_Fc) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{0.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{0.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{3.}
    );

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc =
        CreateView<double[number_of_qps][6]>("qp_Fc", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{6., 12., 18., 24., 30., 36.};
    const auto resid_exact =
        Kokkos::View<double[1][number_of_nodes][6], Kokkos::HostSpace>::const_type(
            resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementOneNodeOneQP_Fd) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{3.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{4.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{0.}
    );

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc = Kokkos::View<double[number_of_qps][6]>("qp_Fc");
    const auto qp_Fd =
        CreateView<double[number_of_qps][6]>("qp_Fd", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{24., 48., 72., 96., 120., 144.};
    const auto resid_exact =
        Kokkos::View<double[1][number_of_nodes][6], Kokkos::HostSpace>::const_type(
            resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementOneNodeOneQP_Fi) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{3.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{4.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{0.}
    );

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc = Kokkos::View<double[number_of_qps][6]>("qp_Fc");
    const auto qp_Fd =
        CreateView<double[number_of_qps][6]>("qp_Fd", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{24., 48., 72., 96., 120., 144.};
    const auto resid_exact =
        Kokkos::View<double[1][number_of_nodes][6], Kokkos::HostSpace>::const_type(
            resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementOneNodeOneQP_Fe) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{3.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{4.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{0.}
    );

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc = Kokkos::View<double[number_of_qps][6]>("qp_Fc");
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe =
        CreateView<double[number_of_qps][6]>("qp_Fe", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{-24., -48., -72., -96., -120., -144.};
    const auto resid_exact =
        Kokkos::View<const double[1][number_of_nodes][6], Kokkos::HostSpace>(resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementOneNodeOneQP_Fg) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{3.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{4.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{0.}
    );

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc = Kokkos::View<double[number_of_qps][6]>("qp_Fc");
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg =
        CreateView<double[number_of_qps][6]>("qp_Fg", std::array{1., 2., 3., 4., 5., 6.});

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{-24., -48., -72., -96., -120., -144.};
    const auto resid_exact =
        Kokkos::View<const double[1][number_of_nodes][6], Kokkos::HostSpace>(resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementOneNodeOneQP_FX) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{3.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{4.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{0.}
    );

    const auto node_FX =
        CreateView<double[number_of_nodes][6]>("node_FX", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fc = Kokkos::View<double[number_of_qps][6]>("qp_Fc");
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{-1., -2., -3., -4., -5., -6.};
    const auto resid_exact =
        Kokkos::View<const double[1][number_of_nodes][6], Kokkos::HostSpace>(resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, TwoElementsOneNodeOneQP) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{0.});
    const auto shape_interp = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_interp", std::array<double, 1>{0.}
    );
    const auto shape_deriv = CreateLeftView<double[number_of_nodes][number_of_qps]>(
        "shape_deriv", std::array<double, 1>{3.}
    );

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc_1 =
        CreateView<double[number_of_qps][6]>("qp_Fc_1", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fc_2 =
        CreateView<double[number_of_qps][6]>("qp_Fc_2", std::array{2., 4., 6., 8., 10., 12.});
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[2][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_deriv,
            .node_FX_ = node_FX,
            .qp_Fc_ = qp_Fc_1,
            .qp_Fd_ = qp_Fd,
            .qp_Fi_ = qp_Fi,
            .qp_Fe_ = qp_Fe,
            .qp_Fg_ = qp_Fg,
            .residual_vector_terms_ = residual_vector_terms
        }
    );

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 1U,
            .num_qps = number_of_qps,
            .qp_weight_ = qp_weights,
            .qp_jacobian_ = qp_jacobian,
            .shape_interp_ = shape_interp,
            .shape_deriv_ = shape_deriv,
            .node_FX_ = node_FX,
            .qp_Fc_ = qp_Fc_2,
            .qp_Fd_ = qp_Fd,
            .qp_Fi_ = qp_Fi,
            .qp_Fe_ = qp_Fe,
            .qp_Fg_ = qp_Fg,
            .residual_vector_terms_ = residual_vector_terms
        }
    );

    constexpr auto resid_exact_data =
        std::array<double, 2 * number_of_nodes * 6>{6.,  12., 18., 24., 30., 36.,
                                                    12., 24., 36., 48., 60., 72.};
    const auto resid_exact =
        Kokkos::View<const double[2][number_of_nodes][6], Kokkos::HostSpace>(resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementOneNodeTwoQPs) {
    constexpr auto number_of_nodes = size_t{1U};
    constexpr auto number_of_qps = size_t{2U};

    const auto qp_weights = CreateView<double[number_of_qps]>("qp_weights", std::array{2., 1.});
    const auto qp_jacobian = CreateView<double[number_of_qps]>("qp_jacobian", std::array{0., 0.});
    const auto shape_interp =
        CreateLeftView<double[number_of_nodes][number_of_qps]>("shape_interp", std::array{0., 0.});
    const auto shape_deriv =
        CreateLeftView<double[number_of_nodes][number_of_qps]>("shape_deriv", std::array{3., 2.});

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    // Note: using std::vector because of compiler bug in nvcc
    const auto qp_Fc = CreateView<double[number_of_qps][6]>(
        "qp_FC", std::vector{1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.}
    );
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");

    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, number_of_nodes * 6>{20., 28., 36., 44., 52., 60.};
    const auto resid_exact =
        Kokkos::View<const double[1][number_of_nodes][6], Kokkos::HostSpace>(resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

TEST(IntegrateResidualVector, OneElementTwoNodesOneQP) {
    constexpr auto number_of_nodes = size_t{2U};
    constexpr auto number_of_qps = size_t{1U};

    const auto qp_weights =
        CreateView<double[number_of_qps]>("qp_weights", std::array<double, 1>{2.});
    const auto qp_jacobian =
        CreateView<double[number_of_qps]>("qp_jacobian", std::array<double, 1>{0.});
    const auto shape_interp =
        CreateLeftView<double[number_of_nodes][number_of_qps]>("shape_interp", std::array{0., 0.});
    const auto shape_deriv =
        CreateLeftView<double[number_of_nodes][number_of_qps]>("shape_deriv", std::array{3., 0.});

    const auto node_FX = Kokkos::View<double[number_of_nodes][6]>("node_FX");
    const auto qp_Fc =
        CreateView<double[number_of_qps][6]>("qp_FC", std::array{1., 2., 3., 4., 5., 6.});
    const auto qp_Fd = Kokkos::View<double[number_of_qps][6]>("qp_Fd");
    const auto qp_Fi = Kokkos::View<double[number_of_qps][6]>("qp_Fi");
    const auto qp_Fe = Kokkos::View<double[number_of_qps][6]>("qp_Fe");
    const auto qp_Fg = Kokkos::View<double[number_of_qps][6]>("qp_Fg");

    const auto residual_vector_terms =
        Kokkos::View<double[1][number_of_nodes][6]>("residual_vector_terms");
    Kokkos::parallel_for(
        "IntegrateResidualVectorElement", number_of_nodes,
        beams::IntegrateResidualVectorElement<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_qps = number_of_qps,
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

    constexpr auto resid_exact_data =
        std::array<double, 2 * number_of_nodes * 6>{6., 12., 18., 24., 30., 36.,
                                                    0., 0.,  0.,  0.,  0.,  0.};
    const auto resid_exact =
        Kokkos::View<const double[1][number_of_nodes][6], Kokkos::HostSpace>(resid_exact_data.data()
        );

    const auto residual_vector_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_vector_terms);
    CompareWithExpected(residual_vector_terms_mirror, resid_exact);
}

}  // namespace kynema::beams::tests
