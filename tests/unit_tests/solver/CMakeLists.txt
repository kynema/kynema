# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_compute_col_inds.cpp
    test_compute_row_ptrs.cpp
)
