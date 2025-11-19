# Specify the source files for the regression test executable
target_sources(
    kynema_regression_tests
    PRIVATE
    dynamic_composite_beam_bending.cpp
    dynamic_IEA15MW_blade_bending.cpp
    dynamic_heavy_top.cpp
    static_composite_beam_bending.cpp
    static_curved_composite_beam_bending.cpp
    static_twisted_composite_beam_bending.cpp
    static_isotropic_beam_rollup.cpp
)
