# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_netcdf_file.cpp
    test_node_state_writer.cpp
    test_time_series_writer.cpp
)
