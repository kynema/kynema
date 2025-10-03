#include <array>
#include <ranges>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "math/vector_operations.hpp"

namespace kynema::tests {

template <unsigned size>
Kokkos::View<double[size]> Create1DView(const std::array<double, size>& input) {
    auto view = Kokkos::View<double[size]>("view");
    auto view_host =
        typename Kokkos::View<double[size], Kokkos::HostSpace>::const_type(input.data());
    Kokkos::deep_copy(view, view_host);
    return view;
}

Kokkos::View<double[3][3]> TestVecTilde(const Kokkos::View<double[3]>& v) {
    auto m = Kokkos::View<double[3][3]>("m");
    Kokkos::parallel_for("VecTilde", 1, KOKKOS_LAMBDA(int) { math::VecTilde(v, m); });
    return m;
}

TEST(VectorTest, VecTilde) {
    auto v = Create1DView<3>({1., 2., 3.});
    const auto m = TestVecTilde(v);

    const auto m_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m);

    constexpr auto expected_data = std::array{0., -3., 2., 3., 0., -1., -2., 1., 0.};
    const auto expected =
        Kokkos::View<double[3][3], Kokkos::HostSpace>::const_type(expected_data.data());

    for (auto i : std::views::iota(0, 3)) {
        for (auto j : std::views::iota(0, 3)) {
            EXPECT_NEAR(m_mirror(i, j), expected(i, j), 1.e-15);
        }
    }
}

void test_DotProduct_View() {
    auto a = Create1DView<3>({1., 2., 3.});
    auto b = Create1DView<3>({4., 5., 6.});
    auto c = 0.;
    Kokkos::parallel_reduce(
        "DotProduct_View", 1,
        KOKKOS_LAMBDA(int, double& result) { result = math::DotProduct(a, b); }, c
    );
    ASSERT_EQ(c, 32.);
}

TEST(VectorTest, DotProduct_View) {
    test_DotProduct_View();
}

}  // namespace kynema::tests
