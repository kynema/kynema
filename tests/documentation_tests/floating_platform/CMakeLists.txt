cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project declaration
project(
  FloatingPlatform
  VERSION 0.0.0
  DESCRIPTION "A flexible multibody structural dynamics code for wind turbines"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Kynema REQUIRED)

add_executable(floating_platform)
target_link_libraries(floating_platform PRIVATE KynemaFMB::kynema_fmb_library)
target_sources(floating_platform PRIVATE floating_platform.cpp)
