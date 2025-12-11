cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard to be used")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project declaration
project(
  NREL5MW_Turbine
  VERSION 1.0.0
  DESCRIPTION "Assembles and solves a full turbine assembly with defined torque and pitch control"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Kynema REQUIRED)

add_executable(nrel5mw_turbine)
target_link_libraries(nrel5mw_turbine PRIVATE Kynema::kynema_library)
target_sources(nrel5mw_turbine PRIVATE nrel5mw_turbine.cpp)
