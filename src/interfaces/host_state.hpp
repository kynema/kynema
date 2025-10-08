#pragma once

#include "state/state.hpp"

namespace kynema::interfaces {

/**
 * @brief Host-side mirror of the simulation state for a given time increment
 *
 * @details This struct maintains host-local copies of the key state variables required
 * for the dynamic simulation, including:
 * - Position: Current position and orientation [x, y, z, qw, qx, qy, qz] -> 7 x 1
 * - Displacement: Current displacement and rotation [δx, δy, δz, δqw, δqx, δqy, δqz] -> 7 x 1
 * - Velocity: Current linear and angular velocities [vx, vy, vz, wx, wy, wz] -> 6 x 1
 * - Acceleration: Current linear and angular accelerations [ax, ay, az, αx, αy, αz] -> 6 x 1
 *
 * @note This struct serves as a host-side mirror of the device state, allowing for
 * efficient data transfer between device and host memory. It's primarily used for
 * the interfaces, I/O operations etc.
 */
template <typename DeviceType>
struct HostState {
    template <typename ValueType>
    using HostView = typename Kokkos::View<ValueType, DeviceType>::HostMirror;

    /// @brief Host local copy of current position
    HostView<double* [7]> x;

    /// @brief Host local copy of current displacement
    HostView<double* [7]> q;

    /// @brief Host local copy of current velocity
    HostView<double* [6]> v;

    /// @brief Host local copy of current acceleration
    HostView<double* [6]> vd;

    /// @brief Host local copy of external forces
    HostView<double* [6]> f;

    /// @brief  Construct host state from state
    explicit HostState(const State<DeviceType>& state)
        : x(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.x)),
          q(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.q)),
          v(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.v)),
          vd(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.vd)),
          f(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, state.f)) {}

    /// @brief Copy state data to host state
    void CopyFromState(const State<DeviceType>& state) const {
        Kokkos::deep_copy(this->x, state.x);
        Kokkos::deep_copy(this->q, state.q);
        Kokkos::deep_copy(this->v, state.v);
        Kokkos::deep_copy(this->vd, state.vd);
    }

    /// @brief Set forces in state from host state
    /// @param state State to update
    void CopyForcesToState(State<DeviceType>& state) const { Kokkos::deep_copy(state.f, this->f); }
};

}  // namespace kynema::interfaces
