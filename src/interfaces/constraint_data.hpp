#pragma once

#include <array>
#include <cstddef>

#include <Kokkos_Core.hpp>

namespace kynema::interfaces {

template <typename DeviceType>
struct HostConstraints;

/**
 * @brief A wrapper around the Constraint ID of a given constraint
 */
struct ConstraintData {
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /// @brief Constraint identifier in model
    size_t id;

    /// @brief Constraint inputs
    std::array<double, 6> inputs{0., 0., 0., 0., 0., 0.};

    /// @brief Constraint outputs
    std::array<double, 3> outputs{0., 0., 0.};

    /// @brief Point loads/moment applied by constraint in global coordinates
    std::array<double, 6> loads{0., 0., 0., 0., 0., 0.};

    /// @brief Constraint data constructor
    /// @param id_ Constraint identifier in model
    explicit ConstraintData(size_t id_) : id(id_) {}

    /// @brief Populates constraint Loads from state data
    /// @param host_constraints Host constraints from which to obtain data
    void GetLoads(const HostConstraints<DeviceType>& host_constraints);

    /// @brief Populates constraint outputs from state data
    /// @param host_constraints Host constraints from which to obtain data
    void GetOutputs(const HostConstraints<DeviceType>& host_constraints);

    /// @brief Sets the constraint inputs in the host constraints
    /// @param host_constraints Host constraints to update
    void SetInputs(HostConstraints<DeviceType>& host_constraints) const;
};

}  // namespace kynema::interfaces
