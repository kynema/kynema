configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

add_subdirectory(regression_tests)
add_subdirectory(unit_tests)
