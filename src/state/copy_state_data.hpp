#pragma once

#include "state.hpp"

namespace kynema_fmb {

/**
 * @brief Performs a deep copy of the state data which might have changed in a given
 * time step.
 *
 * @details It is assumed that the target State object has all of its data which
 * is unchanged (connectivity information, x0, ID) already copied over and all of its
 * Views are properly sized.  One way to ensure this is to first create it with the
 * Clone state method.
 *
 * @tparam DeviceType The Kokkos Device where copy and old states reside
 * @param old The State from which to be copied
 * @param copy The State to which to be copied
 */
template <typename DeviceType>
inline void CopyStateData(State<DeviceType>& copy, const State<DeviceType>& old) {
    using Kokkos::deep_copy;

    copy.time_step = old.time_step;
    deep_copy(copy.x, old.x);
    deep_copy(copy.q_delta, old.q_delta);
    deep_copy(copy.q_prev, old.q_prev);
    deep_copy(copy.q, old.q);
    deep_copy(copy.v, old.v);
    deep_copy(copy.vd, old.vd);
    deep_copy(copy.a, old.a);
    deep_copy(copy.f, old.f);
    deep_copy(copy.tangent, old.tangent);
}

}  // namespace kynema_fmb
