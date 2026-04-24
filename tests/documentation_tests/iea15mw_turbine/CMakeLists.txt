cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project declaration
project(
  IEA15MW_Trubine
  VERSION 1.0.0
  DESCRIPTION "Assembles and solves a full turbine assembly with defined torque and pitch control"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Kynema REQUIRED)

add_executable(iea15mw_turbine)
target_link_libraries(iea15mw_turbine PRIVATE KynemaFMB::kynema_fmb_library)
target_sources(iea15mw_turbine PRIVATE iea15mw_turbine.cpp)
