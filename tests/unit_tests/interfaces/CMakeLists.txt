add_subdirectory(components)

# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_outputs.cpp
)
