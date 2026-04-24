#pragma once

#include <istream>

#include <Kokkos_Core.hpp>

#include "state.hpp"

namespace kynema_fmb {

/**
 * @brief Reads State data from a provided restart file.
 *
 * @details It is assumed that the State data structure has a properly initialized
 * connectivity information, such as when it has been created through the Model interface.
 *
 * @tparam DeviceType the Kokkos Device where the State object resides
 *
 * @param input An input stream from a binary restart file
 * @param state The State object into which to write data
 */
template <typename DeviceType>
inline void ReadStateFromFile(std::istream& input, State<DeviceType>& state) {
    auto num_system_nodes = size_t{};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    input.read(reinterpret_cast<char*>(&num_system_nodes), sizeof(size_t));

    if (num_system_nodes != state.num_system_nodes) {
        throw std::length_error("Number of system nodes in file is not the same as in model");
    }

    const auto mirror_7 = Kokkos::View<double* [7]>::HostMirror("mirror_7", num_system_nodes);
    const auto out_7 = Kokkos::View<double* [7], Kokkos::HostSpace>("out_7", num_system_nodes);

    const auto mirror_6 = Kokkos::View<double* [6]>::HostMirror("mirror_6", num_system_nodes);
    const auto out_6 = Kokkos::View<double* [6], Kokkos::HostSpace>("out_6", num_system_nodes);

    const auto read_7 = [&](const Kokkos::View<double* [7]>& data) {
        const auto stream_size = static_cast<long>(7U * num_system_nodes * sizeof(double));
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        input.read(reinterpret_cast<char*>(out_7.data()), stream_size);

        Kokkos::deep_copy(mirror_7, out_7);
        Kokkos::deep_copy(data, mirror_7);
    };

    const auto read_6 = [&](const Kokkos::View<double* [6]>& data) {
        const auto stream_size = static_cast<long>(6U * num_system_nodes * sizeof(double));
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        input.read(reinterpret_cast<char*>(out_6.data()), stream_size);

        Kokkos::deep_copy(mirror_6, out_6);
        Kokkos::deep_copy(data, mirror_6);
    };

    read_7(state.x0);
    read_7(state.x);
    read_6(state.q_delta);
    read_7(state.q_prev);
    read_7(state.q);
    read_6(state.v);
    read_6(state.vd);
    read_6(state.a);
    read_6(state.f);
}

}  // namespace kynema_fmb
