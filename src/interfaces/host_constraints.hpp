#pragma once

#include "constraints/constraints.hpp"

namespace kynema::interfaces {

/**
 * @brief Host-side mirror of the constraint input, output, and loads for a given time increment
 *
 * @details This struct maintains host-local copies of the key constraint values used
 * for the dynamic simulation, including:
 * - Input: Control inputs for the constraints [x, y, z, qw, qx, qy, qz] -> 7 x 1
 * - Output: Output values depending on constraint [δx, δy, δz, δqw, δqx, δqy, δqz] -> 7 x 1
 * - Loads: Constraint loads (Lambda) [fx, fy, fz, mx, my, mz] -> 6 x 1
 *
 * @note This struct serves as a host-side mirror of the constraint data, allowing for
 * efficient data transfer between device and host memory. It's primarily used for
 * the interfaces, I/O operations etc.
 */
template <typename DeviceType>
struct HostConstraints {
    template <typename ValueType>
    using HostView = typename Kokkos::View<ValueType, DeviceType>::HostMirror;

    /// @brief Host local copy of current inputs
    HostView<double* [7]> input;

    /// @brief Host local copy of current outputs
    HostView<double* [3]> output;

    /// @brief Host local copy of current loads
    HostView<double* [6]> loads;

    /// @brief  Construct host state from state
    explicit HostConstraints(const Constraints<DeviceType>& constraints)
        : input(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, constraints.input)),
          output(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, constraints.output)),
          loads(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, constraints.lambda)) {}

    /// @brief Copy state data to host state
    void CopyFromConstraints(const Constraints<DeviceType>& constraints) {
        Kokkos::deep_copy(this->input, constraints.input);
        Kokkos::deep_copy(this->output, constraints.output);
        Kokkos::deep_copy(this->loads, constraints.lambda);
    }

    /// @brief Set inputs in constraints from host constraints
    /// @param constraints constraints to update
    void CopyInputsToConstraints(Constraints<DeviceType>& constraints) const {
        Kokkos::deep_copy(constraints.input, this->input);
    }
};

}  // namespace kynema::interfaces
