# Specify the source files for the regression test executable
target_sources(
    kynema_regression_tests
    PRIVATE
    test_beams.cpp
    test_cantilever_beam.cpp
    test_controller.cpp
    test_heavy_top.cpp
    test_masses.cpp
    test_precession.cpp
    test_rotating_beam.cpp
    test_rotor.cpp
    test_solver.cpp
    test_springs.cpp
    test_spring_mass_system.cpp
    test_utilities.cpp
    test_netcdf_outputs_writer.cpp
)

add_subdirectory(verification)
