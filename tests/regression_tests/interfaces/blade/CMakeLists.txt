# Specify the source files for the regression test executable
target_sources(
    kynema_regression_tests
    PRIVATE
    test_blade_interface.cpp
    test_static_curled_beam_trapezoidal_quadrature.cpp
    test_static_curled_beam_gll_quadrature.cpp
    test_static_curled_beam_gl_quadrature.cpp
)
