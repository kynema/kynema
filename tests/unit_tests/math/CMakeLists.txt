# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_interpolation.cpp
    test_least_squares_fit.cpp
    test_matrix_operations.cpp
    test_project_points_to_target_polynomial.cpp
    test_quaternion_operations.cpp
    test_vector_operations.cpp
)
