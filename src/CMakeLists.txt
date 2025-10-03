#--------------------------------------------------------------------------
# Define the main library target
#--------------------------------------------------------------------------
add_library(kynema_library)
add_library(Kynema::kynema_library ALIAS kynema_library)

#----------------------------------------
# Set target properties
#----------------------------------------
set_target_properties(kynema_library PROPERTIES
  LINKER_LANGUAGE CXX
  POSITION_INDEPENDENT_CODE ON  # Ensure position-independent code for shared libraries
  VISIBILITY_INLINES_HIDDEN YES # Hide inline function symbols to reduce binary size
)

#----------------------------------------
# Link internal libraries and dependencies
#----------------------------------------
target_link_libraries(kynema_library PRIVATE
  kynema_options
  kynema_warnings
)

#----------------------------------------
# Add compile definitions
#----------------------------------------
if(Kynema_WRITE_OUTPUTS)
  target_compile_definitions(kynema_library PRIVATE Kynema_WRITE_OUTPUTS)
endif()

#--------------------------------------------------------------------------
# Link system libraries
#--------------------------------------------------------------------------
#----------------------------------------
# yaml-cpp
#----------------------------------------
if(TARGET yaml-cpp::yaml-cpp)
  set(YAML_CPP_TARGET yaml-cpp::yaml-cpp)
else()
  set(YAML_CPP_TARGET yaml-cpp)
endif()

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

#--------------------------------------------------------------------------
# Include directories for system libraries
#--------------------------------------------------------------------------
target_link_system_libraries(kynema_library PUBLIC
  KokkosKernels::kokkoskernels
  LAPACK::LAPACK
  Eigen3::Eigen
  ${YAML_CPP_TARGET}
  ${FS_LIB}
  ${NETCDF_TARGET}
)

if(Kynema_ENABLE_CUDSS)
target_link_system_libraries(kynema_library PUBLIC
  cudss
)
endif()

if(Kynema_ENABLE_MKL)
target_link_system_libraries(kynema_library PUBLIC
  MKL::MKL
)
endif()

if(Kynema_ENABLE_SUPERLU)
target_link_system_libraries(kynema_library PUBLIC
  superlu::superlu
)
endif()

if(Kynema_ENABLE_SUPERLU_MT)
target_link_system_libraries(kynema_library PUBLIC
  superlu_mt::superlu_mt
)
endif()

if(Kynema_ENABLE_KLU)
target_link_system_libraries(kynema_library PUBLIC
  SuiteSparse::KLU
)
endif()

if(Kynema_ENABLE_UMFPACK)
target_link_system_libraries(kynema_library PUBLIC
  SuiteSparse::UMFPACK
)
endif()

#--------------------------------------------------------------------------
# Include shared libraries
#--------------------------------------------------------------------------
include(../cmake/SharedLibraries.cmake)

#----------------------------------------
# Define shared libraries for controllers
#----------------------------------------
add_shared_library(DISCON)
add_shared_library(DISCON_ROTOR_TEST_CONTROLLER)

#----------------------------------------
# Set include directories
#----------------------------------------
if(Kynema_ENABLE_MKL)
else()
find_path(LAPACKE_INCLUDE_DIRS lapacke.h)
target_include_directories(kynema_library SYSTEM
 PUBLIC
${LAPACKE_INCLUDE_DIRS}
)
endif()

target_include_directories(kynema_library
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/Kynema>
)

# Set the C++ standard for the target
target_compile_features(kynema_library PUBLIC cxx_std_20)

#----------------------------------------
# Add subdirectories
#----------------------------------------
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
add_subdirectory(vendor)

#----------------------------------------
# Install targets and export configuration
#----------------------------------------
install(
  TARGETS kynema_library kynema_options kynema_warnings
  EXPORT KynemaTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${PACKAGE_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT KynemaTargets
  FILE KynemaTargets.cmake
  NAMESPACE Kynema::
  DESTINATION lib/cmake/Kynema
)

install(FILES
  DESTINATION include/Kynema/
)
