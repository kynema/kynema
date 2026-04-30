# Test ROSCO controller
if (KYNEMA_FMB_ENABLE_ROSCO_CONTROLLER AND
    (NOT KYNEMA_FMB_ENABLE_SANITIZER_LEAK))
    target_sources(
        kynema_fmb_regression_tests
        PRIVATE
        test_rosco_controller.cpp
        test_turbine_rosco_interface.cpp
    )
endif()

# Test AeroDyn Inflow
if (KYNEMA_FMB_ENABLE_OPENFAST_ADI AND
    KYNEMA_FMB_ENABLE_ROSCO_CONTROLLER AND
    (NOT KYNEMA_FMB_ENABLE_SANITIZER_LEAK))
    target_sources(
        kynema_fmb_regression_tests
        PRIVATE
        test_rotor_aero_controller.cpp
    )
endif()

add_custom_target(
    copy_external_test_files ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/IEA-15-240-RWT
            $<TARGET_FILE_DIR:kynema_fmb_regression_tests>/IEA-15-240-RWT)
