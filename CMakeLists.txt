# -----------------------------------------------------------------------------
# Basic project setup
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

# Set C++ standard requirements
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project declaration
project(
  Kynema
  VERSION 0.9.0
  DESCRIPTION "A performance-oriented flexible multibody structural dynamics code"
  LANGUAGES CXX
)

# -----------------------------------------------------------------------------
# Project configuration and options
# -----------------------------------------------------------------------------
# Include custom CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(cmake/PreventInSourceBuilds.cmake)
include(cmake/KynemaOptions.cmake)

# Setup project-specific options
kynema_setup_options()
kynema_global_options()

# Include and setup dependencies
include(cmake/Dependencies.cmake)
kynema_setup_dependencies()

# Setup local project options
kynema_local_options()

# -----------------------------------------------------------------------------
# Version information
# -----------------------------------------------------------------------------
# Generate a shortened Git SHA for versioning
set(GIT_SHA "Unknown" CACHE STRING "SHA this build was generated from")
string(SUBSTRING "${GIT_SHA}" 0 8 GIT_SHORT_SHA)

# Set Kynema Version
set(KYNEMA_VERSION_MAJOR 0 CACHE STRING "major version" FORCE)
set(KYNEMA_VERSION_MINOR 1 CACHE STRING "minor version" FORCE)
set(KYNEMA_VERSION ${KYNEMA_VERSION_MAJOR}.${KYNEMA_VERSION_MINOR} CACHE STRING "version" FORCE)

# -----------------------------------------------------------------------------
# Library targets
# -----------------------------------------------------------------------------
# Alias libraries for ease of use in targets
add_library(Kynema::kynema_options ALIAS kynema_options)
add_library(Kynema::kynema_warnings ALIAS kynema_warnings)

# Add the main source directory
add_subdirectory(src)

# -----------------------------------------------------------------------------
# Installation configuration
# -----------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/KynemaConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/KynemaConfig.cmake"
  INSTALL_DESTINATION lib/cmake/Kynema
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/KynemaConfigVersion.cmake"
  VERSION "${version}"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/KynemaConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/KynemaConfigVersion.cmake"
  DESTINATION lib/cmake/Kynema
)

# -----------------------------------------------------------------------------
# Optional components
# -----------------------------------------------------------------------------
# Build tests
include(CTest)
if(Kynema_ENABLE_TESTS)
    add_subdirectory(tests)
endif()

# Build documentation
if(Kynema_ENABLE_DOCUMENTATION)
    add_subdirectory(docs)
endif()

configure_file(cmake/Kynema_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/Kynema_config.h)
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/src/Kynema_config.h
        DESTINATION include/Kynema/
)
