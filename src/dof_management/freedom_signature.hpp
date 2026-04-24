#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

namespace kynema_fmb::dof {

/**
 * @brief Represents the active degrees of freedom for a node
 *
 * Each bit in the 6-bit signature represents a degree of freedom:
 * - Bits 5-3: Translational DOFs (x, y, z positions)
 * - Bits 2-0: Rotational DOFs (rx, ry, rz rotations)
 *
 * Binary format: 0b00[Tz][Ty][Tx][Rz][Ry][Rx]
 */
enum class FreedomSignature : std::uint8_t {
    AllComponents = 0b00111111,  //< Enable all 6 degrees of freedom
    JustPosition = 0b00111000,   //< Enable only translational DOFs
    JustRotation = 0b00000111,   //< Enable only rotational DOFs
    NoComponents = 0b00000000    //< Disable all degrees of freedom
};

/**
 * @brief Combines two freedom signatures using bitwise OR
 *
 * @param x First freedom signature
 * @param y Second freedom signature
 * @return Combined freedom signature with all DOFs from both inputs enabled
 */
KOKKOS_INLINE_FUNCTION
FreedomSignature operator|(FreedomSignature x, FreedomSignature y) {
    using T = std::underlying_type_t<FreedomSignature>;
    return static_cast<FreedomSignature>(static_cast<T>(x) | static_cast<T>(y));
}

/**
 * @brief Counts the number of active degrees of freedom in a signature
 *
 * @param x The freedom signature to count
 * @return The number of active degrees of freedom
 */
KOKKOS_INLINE_FUNCTION
size_t count_active_dofs(FreedomSignature x) {
    using T = std::underlying_type_t<FreedomSignature>;
    auto count = 0UL;
    constexpr auto zero = T{0};
    constexpr auto one = T{1};
    for (auto value = static_cast<T>(x); value > zero; value = value >> 1) {
        count += value & one;
    }
    return count;
}

}  // namespace kynema_fmb::dof
