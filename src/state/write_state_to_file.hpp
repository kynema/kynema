#pragma once

#include <ostream>

#include <Kokkos_Core.hpp>

#include "state.hpp"

namespace kynema_fmb {

/**
 * @brief Writes State data into a minimal restart file
 *
 * @tparam DeviceType the Kokkos Device where the State object resides
 *
 * @param output An output stream to the binary restart file
 * @param state The State object from which to write out the data
 */
template <typename DeviceType>
inline void WriteStateToFile(std::ostream& output, const State<DeviceType>& state) {
    auto num_system_nodes = state.num_system_nodes;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    output.write(reinterpret_cast<char*>(&num_system_nodes), sizeof(size_t));

    const auto mirror_7 = Kokkos::View<double* [7]>::HostMirror("mirror_7", num_system_nodes);
    const auto out_7 = Kokkos::View<double* [7], Kokkos::HostSpace>("out_7", num_system_nodes);

    const auto mirror_6 = Kokkos::View<double* [6]>::HostMirror("mirror_6", num_system_nodes);
    const auto out_6 = Kokkos::View<double* [6], Kokkos::HostSpace>("out_6", num_system_nodes);

    const auto write_7 = [&](const Kokkos::View<double* [7]>& data) {
        Kokkos::deep_copy(mirror_7, data);
        Kokkos::deep_copy(out_7, mirror_7);

        const auto stream_size = static_cast<long>(7U * num_system_nodes * sizeof(double));
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        output.write(reinterpret_cast<char*>(out_7.data()), stream_size);
    };

    const auto write_6 = [&](const Kokkos::View<double* [6]>& data) {
        Kokkos::deep_copy(mirror_6, data);
        Kokkos::deep_copy(out_6, mirror_6);

        const auto stream_size = static_cast<long>(6U * num_system_nodes * sizeof(double));
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        output.write(reinterpret_cast<char*>(out_6.data()), stream_size);
    };

    write_7(state.x0);
    write_7(state.x);
    write_6(state.q_delta);
    write_7(state.q_prev);
    write_7(state.q);
    write_6(state.v);
    write_6(state.vd);
    write_6(state.a);
    write_6(state.f);
}

}  // namespace kynema_fmb
