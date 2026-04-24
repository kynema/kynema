cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project declaration
project(
  ThreeBladeRotor
  VERSION 1.0.0
  DESCRIPTION ""
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Kynema REQUIRED)

add_executable(three_blade_rotor)
target_link_libraries(three_blade_rotor PRIVATE KynemaFMB::kynema_fmb_library)
target_sources(three_blade_rotor PRIVATE three_blade_rotor.cpp)
