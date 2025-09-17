#include <numbers>
#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "math/matrix_operations.hpp"

namespace kynema::tests {

template <size_t rows, size_t cols>
Kokkos::View<double[rows][cols]> Create2DView(const std::array<double, rows * cols>& input) {
    auto view =
        Kokkos::View<double[rows][cols]>(Kokkos::view_alloc("view", Kokkos::WithoutInitializing));
    auto view_host = Kokkos::View<const double[rows][cols], Kokkos::HostSpace>(input.data());
    auto view_mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, view);
    Kokkos::deep_copy(view_mirror, view_host);
    Kokkos::deep_copy(view, view_mirror);
    return view;
}

inline void test_AX_Matrix() {
    const auto A = Create2DView<3, 3>({1., 2., 3., 4., 5., 6., 7., 8., 9.});
    const auto out = Kokkos::View<double[3][3]>("out");
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) { math::AX_Matrix(A, out); });
    const auto out_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);

    constexpr auto expected_data = std::array{7., -1., -1.5, -2., 5., -3., -3.5, -4., 3.};
    const auto expected =
        Kokkos::View<double[3][3], Kokkos::HostSpace>::const_type(expected_data.data());

    for (auto i : std::views::iota(0, 3)) {
        for (auto j : std::views::iota(0, 3)) {
            EXPECT_NEAR(out_mirror(i, j), expected(i, j), 1.e-15);
        }
    }
}

TEST(MatrixTest, AX_Matrix) {
    test_AX_Matrix();
}

Kokkos::View<double[3]> TestAxialVectorOfMatrix(const Kokkos::View<const double[3][3]>& m) {
    auto v = Kokkos::View<double[3]>("v");
    Kokkos::parallel_for(
        "AxialVectorOfMatrix", 1, KOKKOS_LAMBDA(int) { math::AxialVectorOfMatrix(m, v); }
    );
    return v;
}

TEST(MatrixTest, AxialVectorOfMatrix) {
    const auto m = Create2DView<3, 3>({0., -1., 0., 1., 0., 0., 0., 0., 0.});
    const auto v = TestAxialVectorOfMatrix(m);

    const auto v_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);

    constexpr auto expected_data = std::array{0., 0., 1.};
    const auto expected =
        Kokkos::View<double[3], Kokkos::HostSpace>::const_type(expected_data.data());

    for (auto i : std::views::iota(0, 3)) {
        EXPECT_NEAR(v_mirror(i), expected(i), 1.e-15);
    }
}

TEST(MatrixTest, RotateMatrix6_1) {
    const auto m = std::array{
        std::array{1., 2., 3., 4., 5., 6.},       std::array{7., 8., 9., 10., 11., 12.},
        std::array{13., 14., 15., 16., 17., 18.}, std::array{19., 20., 21., 22., 23., 24.},
        std::array{25., 26., 27., 28., 29., 30.}, std::array{31., 32., 33., 34., 35., 36.},
    };
    const auto m_exp = std::array{
        std::array{1., 2., 3., 4., 5., 6.},       std::array{7., 8., 9., 10., 11., 12.},
        std::array{13., 14., 15., 16., 17., 18.}, std::array{19., 20., 21., 22., 23., 24.},
        std::array{25., 26., 27., 28., 29., 30.}, std::array{31., 32., 33., 34., 35., 36.},
    };
    const auto q = std::array{1., 0., 0., 0.};
    const auto m_act = math::RotateMatrix6(m, q);
    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(m_act[i][j], m_exp[i][j], 1.e-12);
        }
    }
}

TEST(MatrixTest, RotateMatrix6_2) {
    const auto m = std::array{
        std::array{1., 2., 3., 4., 5., 6.},       std::array{7., 8., 9., 10., 11., 12.},
        std::array{13., 14., 15., 16., 17., 18.}, std::array{19., 20., 21., 22., 23., 24.},
        std::array{25., 26., 27., 28., 29., 30.}, std::array{31., 32., 33., 34., 35., 36.},
    };
    const auto m_exp = std::array{
        std::array{1, -2, -3, 4, -5, -6},     std::array{-7, 8, 9, -10, 11, 12},
        std::array{-13, 14, 15, -16, 17, 18}, std::array{19, -20, -21, 22, -23, -24},
        std::array{-25, 26, 27, -28, 29, 30}, std::array{-31, 32, 33, -34, 35, 36},
    };
    const auto q = std::array{0., 1., 0., 0.};
    const auto m_act = math::RotateMatrix6(m, q);
    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(m_act[i][j], m_exp[i][j], 1.e-12);
        }
    }
}

TEST(MatrixTest, RotateMatrix6_3) {
    const auto m = std::array{
        std::array{1., 2., 3., 4., 5., 6.},       std::array{7., 8., 9., 10., 11., 12.},
        std::array{13., 14., 15., 16., 17., 18.}, std::array{19., 20., 21., 22., 23., 24.},
        std::array{25., 26., 27., 28., 29., 30.}, std::array{31., 32., 33., 34., 35., 36.},
    };
    const auto m_exp = std::array{
        std::array{
            16., 11.313708498984761, 1.9999999999999982, 22.000000000000004, 15.556349186104047,
            1.9999999999999987
        },
        std::array{
            11.313708498984759, 8., std::numbers::sqrt2, 15.556349186104045, 11., std::numbers::sqrt2
        },
        std::array{
            11.999999999999998, 8.4852813742385695, -1.0638139115010993E-15, 11.999999999999998,
            8.4852813742385695, -2.3198878784481192E-15
        },
        std::array{
            52.000000000000007, 36.76955262170047, 1.999999999999996, 58., 41.012193308819761,
            1.9999999999999949
        },
        std::array{
            36.76955262170047, 26., std::numbers::sqrt2, 41.012193308819754, 29., std::numbers::sqrt2
        },
        std::array{
            11.999999999999996, 8.4852813742385678, -5.4353103904786874E-16, 11.999999999999995,
            8.4852813742385678, -1.7996050059948887E-15
        },

    };
    const auto q = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(std::numbers::pi / 4., Eigen::Matrix<double, 3, 1>::Unit(1))
    );
    const auto q_array = std::array{q.w(), q.x(), q.y(), q.z()};
    const auto m_act = math::RotateMatrix6(m, q_array);
    for (auto i : std::views::iota(0U, 6U)) {
        for (auto j : std::views::iota(0U, 6U)) {
            EXPECT_NEAR(m_act[i][j], m_exp[i][j], 1.e-12);
        }
    }
}

}  // namespace kynema::tests
