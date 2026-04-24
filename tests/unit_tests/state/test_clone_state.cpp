#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "dof_management/freedom_signature.hpp"
#include "state/clone_state.hpp"
#include "state/state.hpp"

namespace {

template <typename T>
void Compare(const T& field_1, const T& field_2) {
    const auto mirror_1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field_1);
    const auto mirror_2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field_2);

    if constexpr (T::rank() == 1) {
        for (auto i : std::views::iota(0U, field_1.extent(0))) {
            EXPECT_EQ(mirror_1(i), mirror_2(i));
        }
    } else if constexpr (T::rank() == 2) {
        for (auto i : std::views::iota(0U, field_1.extent(0))) {
            for (auto j : std::views::iota(0U, field_1.extent(1))) {
                EXPECT_EQ(mirror_1(i, j), mirror_2(i, j));
            }
        }
    } else if constexpr (T::rank() == 3) {
        for (auto i : std::views::iota(0U, field_1.extent(0))) {
            for (auto j : std::views::iota(0U, field_1.extent(1))) {
                for (auto k : std::views::iota(0U, field_1.extent(2))) {
                    EXPECT_EQ(mirror_1(i, j, k), mirror_2(i, j, k));
                }
            }
        }
    }
}

template <typename DeviceType>
void CompareStates(
    const kynema_fmb::State<DeviceType>& state_1, const kynema_fmb::State<DeviceType>& state_2
) {
    EXPECT_EQ(state_1.num_system_nodes, state_2.num_system_nodes);
    Compare(state_1.ID, state_2.ID);
    Compare(state_1.node_freedom_allocation_table, state_2.node_freedom_allocation_table);
    Compare(state_1.node_freedom_map_table, state_2.node_freedom_map_table);
    Compare(state_1.x0, state_2.x0);
    Compare(state_1.x, state_2.x);
    Compare(state_1.q_delta, state_2.q_delta);
    Compare(state_1.q_prev, state_2.q_prev);
    Compare(state_1.q, state_2.q);
    Compare(state_1.v, state_2.v);
    Compare(state_1.vd, state_2.vd);
    Compare(state_1.a, state_2.a);
    Compare(state_1.tangent, state_2.tangent);
}

template <typename DeviceType>
kynema_fmb::State<DeviceType> CreateTestState() {
    constexpr auto num_system_nodes = 2UL;
    auto state = kynema_fmb::State<DeviceType>(num_system_nodes);
    Kokkos::deep_copy(state.ID, 1UL);
    Kokkos::deep_copy(
        state.node_freedom_allocation_table, kynema_fmb::dof::FreedomSignature::AllComponents
    );
    Kokkos::deep_copy(state.node_freedom_map_table, 3UL);
    Kokkos::deep_copy(state.x0, 4.);
    Kokkos::deep_copy(state.x, 5.);
    Kokkos::deep_copy(state.q_delta, 6.);
    Kokkos::deep_copy(state.q_prev, 7.);
    Kokkos::deep_copy(state.q, 8.);
    Kokkos::deep_copy(state.v, 9.);
    Kokkos::deep_copy(state.vd, 10.);
    Kokkos::deep_copy(state.a, 11.);
    Kokkos::deep_copy(state.tangent, 12.);
    return state;
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CloneState, CloneState) {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
    auto state_1 = CreateTestState<DeviceType>();
    auto state_2 = CloneState(state_1);

    CompareStates(state_1, state_2);
}

}  // namespace kynema_fmb::tests
