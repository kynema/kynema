# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_rotate_section_matrix.cpp
    test_calculate_mass_matrix_components.cpp
    test_calculate_inertial_forces.cpp
    test_calculate_gravity_force.cpp
    test_calculate_gyroscopic_matrix.cpp
    test_calculate_inertia_stiffness_matrix.cpp
)
