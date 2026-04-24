target_sources(kynema_fmb_library PRIVATE)

install(FILES
        calculate_constraint_output.hpp
        calculate_constraint_residual_gradient.hpp
        calculate_fixed_bc_constraint.hpp
        calculate_fixed_bc_3DOF_constraint.hpp
        calculate_prescribed_bc_constraint.hpp
        calculate_prescribed_bc_3DOF_constraint.hpp
        calculate_revolute_joint_constraint.hpp
        calculate_revolute_joint_force.hpp
        calculate_revolute_joint_output.hpp
        calculate_rigid_joint_constraint.hpp
        calculate_rigid_joint_3DOF_constraint.hpp
        calculate_rotation_control_constraint.hpp
        constraint.hpp
        constraints.hpp
        constraint_type.hpp
        update_lambda_prediction.hpp
        DESTINATION include/KynemaFMB/constraints
)
