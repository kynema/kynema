# Specify the source files for the regression test executable
target_sources(
    kynema_regression_tests
    PRIVATE
    static_isotropic_beam_rollup.cpp
    static_composite_beam_bending.cpp
)
