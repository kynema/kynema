#include <gtest/gtest.h>

#include "constraints/constraint_type.hpp"

namespace kynema_fmb::constraints::tests {

TEST(ConstraintTypeTest, NoneConstraintHasOneNode) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::None), 1U);
}

TEST(ConstraintTypeTest, FixedBCHasOneNode) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::FixedBC), 1U);
}

TEST(ConstraintTypeTest, PrescribedBCHasOneNode) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::PrescribedBC), 1U);
}

TEST(ConstraintTypeTest, RigidJointHasTwoNodes) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::RigidJoint), 2U);
}

TEST(ConstraintTypeTest, RevoluteJointHasTwoNodes) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::RevoluteJoint), 2U);
}

TEST(ConstraintTypeTest, RotationControlHasTwoNodes) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::RotationControl), 2U);
}

TEST(ConstraintTypeTest, FixedBC3DOFsHasOneNode) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::FixedBC3DOFs), 1U);
}

TEST(ConstraintTypeTest, PrescribedBC3DOFsHasOneNode) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::PrescribedBC3DOFs), 1U);
}

TEST(ConstraintTypeTest, RigidJoint6DOFsTo3DOFsHasTwoNodes) {
    EXPECT_EQ(GetNumberOfNodes(ConstraintType::RigidJoint6DOFsTo3DOFs), 2U);
}

TEST(ConstraintTypeTest, RevoluteJointHasFiveDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::RevoluteJoint), 5U);
}

TEST(ConstraintTypeTest, NoneConstraintHasSixDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::None), 6U);
}

TEST(ConstraintTypeTest, FixedBCHasSixDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::FixedBC), 6U);
}

TEST(ConstraintTypeTest, PrescribedBCHasSixDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::PrescribedBC), 6U);
}

TEST(ConstraintTypeTest, RigidJointHasSixDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::RigidJoint), 6U);
}

TEST(ConstraintTypeTest, RotationControlHasSixDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::RotationControl), 6U);
}

TEST(ConstraintTypeTest, FixedBC3DOFsHasThreeDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::FixedBC3DOFs), 3U);
}

TEST(ConstraintTypeTest, PrescribedBC3DOFsHasThreeDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::PrescribedBC3DOFs), 3U);
}

TEST(ConstraintTypeTest, RigidJoint6DOFsTo3DOFsHasThreeDOFs) {
    EXPECT_EQ(NumRowsForConstraint(ConstraintType::RigidJoint6DOFsTo3DOFs), 3U);
}

}  // namespace kynema_fmb::constraints::tests
