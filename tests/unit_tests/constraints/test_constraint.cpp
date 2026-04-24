#include <array>
#include <cstddef>

#include <gtest/gtest.h>

#include "constraints/constraint.hpp"
#include "constraints/constraint_type.hpp"

namespace kynema_fmb::constraints::tests {

class ConstraintTest : public ::testing::Test {
protected:
    const size_t id{1};        // constraint id
    const size_t node1_id{1};  // base node id
    const size_t node2_id{2};  // target node id
    const std::array<size_t, 2> node_ids{node1_id, node2_id};
    const std::array<double, 3> axis_vec{0., 1., 0.};  // axis vector is y-axis
    const std::array<double, 7> init_displacement{1.,    2., 3., 0.707,
                                                  0.707, 0., 0.};  // reference displacement
};

TEST_F(ConstraintTest, DefaultInitialization) {
    const auto constraint = Constraint(id, ConstraintType::FixedBC, node_ids);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::FixedBC);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, (std::array{0., 0., 0.}));
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, FixedBCInitialization) {
    const auto constraint = Constraint(id, ConstraintType::FixedBC, node_ids);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::FixedBC);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, (std::array{0., 0., 0.}));
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, RevoluteJointInitialization) {
    const auto constraint = Constraint(id, ConstraintType::RevoluteJoint, node_ids, axis_vec);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::RevoluteJoint);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, PrescribedBCInitialization) {
    const auto constraint =
        Constraint(id, ConstraintType::PrescribedBC, node_ids, axis_vec, init_displacement);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::PrescribedBC);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, init_displacement);
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, RotationControlInitialization) {
    double control_signal = 0.;
    const auto constraint = Constraint(
        id, ConstraintType::RotationControl, node_ids, axis_vec, {0., 0., 0., 1., 0., 0., 0.},
        &control_signal
    );

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::RotationControl);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, &control_signal);
}

TEST_F(ConstraintTest, RigidJointInitialization) {
    const auto constraint = Constraint(id, ConstraintType::RigidJoint, node_ids, axis_vec);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::RigidJoint);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, FixedBC3DOFsInitialization) {
    const auto constraint = Constraint(id, ConstraintType::FixedBC3DOFs, node_ids, axis_vec);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::FixedBC3DOFs);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, PrescribedBC3DOFsInitialization) {
    const auto constraint =
        Constraint(id, ConstraintType::PrescribedBC3DOFs, node_ids, axis_vec, init_displacement);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::PrescribedBC3DOFs);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, init_displacement);
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, RigidJoint6DOFsTo3DOFsInitialization) {
    const auto constraint =
        Constraint(id, ConstraintType::RigidJoint6DOFsTo3DOFs, node_ids, axis_vec);

    EXPECT_EQ(constraint.ID, id);
    EXPECT_EQ(constraint.type, ConstraintType::RigidJoint6DOFsTo3DOFs);
    EXPECT_EQ(constraint.node_ids[0], node1_id);
    EXPECT_EQ(constraint.node_ids[1], node2_id);
    EXPECT_EQ(constraint.axis_vector, axis_vec);
    EXPECT_EQ(constraint.initial_displacement, (std::array{0., 0., 0., 1., 0., 0., 0.}));
    EXPECT_EQ(constraint.control, nullptr);
}

TEST_F(ConstraintTest, CustomInitialDisplacement) {
    const auto custom_disp = std::array{1.5, 2.5, 3.5, 0.866, 0.5, 0., 0.};
    const auto constraint =
        Constraint(id, ConstraintType::PrescribedBC, node_ids, axis_vec, custom_disp);

    EXPECT_EQ(constraint.initial_displacement, custom_disp);
}

TEST_F(ConstraintTest, NullControlSignal) {
    const auto constraint = Constraint(
        id, ConstraintType::RotationControl, node_ids, axis_vec, {0., 0., 0., 1., 0., 0., 0.},
        nullptr
    );

    EXPECT_EQ(constraint.control, nullptr);
}

}  // namespace kynema_fmb::constraints::tests
