# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_calculate_fixed_bc_constraint.cpp
    test_calculate_prescribed_bc_constraint.cpp
    test_calculate_revolute_joint_constraint.cpp
    test_calculate_revolute_joint_force.cpp
    test_calculate_revolute_joint_output.cpp
    test_calculate_rigid_joint_constraint.cpp
    test_calculate_rotation_control_constraint.cpp
    test_constraint.cpp
    test_constraints.cpp
    test_constraint_type.cpp
)
