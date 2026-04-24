cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project declaration
project(
  SpringMassSystem
  VERSION 1.0.0
  DESCRIPTION ""
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Kynema REQUIRED)

add_executable(spring_mass_system)
target_link_libraries(spring_mass_system PRIVATE KynemaFMB::kynema_fmb_library)
target_sources(spring_mass_system PRIVATE spring_mass_system.cpp)
