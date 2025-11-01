#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/masses/calculate_gravity_force.hpp"
#include "test_calculate.hpp"

namespace kynema::tests {

struct ExecuteCalculateGravityForce {
    double mass;
    Kokkos::View<double[3]>::const_type gravity;
    Kokkos::View<double[3][3]>::const_type eta_tilde;
    Kokkos::View<double[6]> FG;

    KOKKOS_FUNCTION
    void operator()(size_t) const {
        masses::CalculateGravityForce<Kokkos::DefaultExecutionSpace>::invoke(
            mass, gravity, eta_tilde, FG
        );
    }
};

TEST(CalculateGravityForceTestsMasses, OneNode) {
    constexpr auto mass = 1.;

    const auto eta_tilde = CreateView<double[3][3]>(
        "eta_tilde", std::array{37., 38., 39., 40., 41., 42., 43., 44., 45.}
    );
    const auto gravity = CreateView<double[3]>("gravity", std::array{46., 47., 48.});

    const auto FG = Kokkos::View<double[6]>("FG");

    Kokkos::parallel_for(
        "CalculateGravityForce", 1,
        ExecuteCalculateGravityForce{
            .mass = mass, .gravity = gravity, .eta_tilde = eta_tilde, .FG = FG
        }
    );

    constexpr auto FG_exact_data = std::array{46., 47., 48., 5360., 5783., 6206.};
    const auto FG_exact = Kokkos::View<const double[6], Kokkos::HostSpace>(FG_exact_data.data());

    const auto FG_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), FG);
    CompareWithExpected(FG_mirror, FG_exact);
}

}  // namespace kynema::tests
