#include <gtest/gtest.h>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::dof::tests {

TEST(TestFreedomSignature, CombineSignatures_NoOverlap) {
    const auto x = FreedomSignature::JustPosition;
    const auto y = FreedomSignature::JustRotation;
    const auto z = x | y;

    EXPECT_EQ(z, FreedomSignature::AllComponents);
}

TEST(TestFreedomSignature, CombineSignatures_Overlap) {
    const auto x = FreedomSignature::JustPosition;
    const auto y = FreedomSignature::AllComponents;
    const auto z = x | y;

    EXPECT_EQ(z, FreedomSignature::AllComponents);
}

TEST(TestFreedomSignature, CountActiveDofs_Position) {
    auto x = count_active_dofs(FreedomSignature::JustPosition);
    EXPECT_EQ(x, 3);
}

TEST(TestFreedomSignature, CountActiveDofs_Rotation) {
    auto x = count_active_dofs(FreedomSignature::JustRotation);
    EXPECT_EQ(x, 3);
}

TEST(TestFreedomSignature, CountActiveDofs_AllComponents) {
    auto x = count_active_dofs(FreedomSignature::AllComponents);
    EXPECT_EQ(x, 6);
}
}  // namespace kynema_fmb::dof::tests
