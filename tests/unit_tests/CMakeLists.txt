# Add the executable target for unit tests
add_executable(kynema_unit_tests)

# Add subdirectories for additional components
add_subdirectory(constraints)
add_subdirectory(dof_management)
add_subdirectory(elements)
add_subdirectory(interfaces)
add_subdirectory(math)
add_subdirectory(model)
add_subdirectory(solver)
add_subdirectory(state)
add_subdirectory(step)
add_subdirectory(system)
add_subdirectory(utilities)

# Specify the source files for the unit tests
target_sources(kynema_unit_tests PRIVATE utest_main.cpp)

target_include_directories(kynema_unit_tests PRIVATE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests/unit_tests>
)

# Link internal libraries to the unit test executable
target_link_libraries(kynema_unit_tests PRIVATE
  kynema_warnings
  kynema_options
  KynemaFMB::kynema_fmb_library
)

#----------------------------------------
# NetCDF
#----------------------------------------
if(TARGET netCDF::netcdf)
  set(NETCDF_TARGET netCDF::netcdf)
elseif(TARGET NetCDF::NetCDF)
  set(NETCDF_TARGET NetCDF::NetCDF)
elseif(TARGET netcdf)
  set(NETCDF_TARGET netcdf)
else()
  message(FATAL_ERROR "NetCDF library not found")
endif()

target_link_system_libraries(kynema_unit_tests PRIVATE
  GTest::gtest
  GTest::gtest_main
)

# Include GoogleTest module and discover tests
include(GoogleTest)
gtest_discover_tests(kynema_unit_tests PROPERTIES DISCOVERY_TIMEOUT 60000)
