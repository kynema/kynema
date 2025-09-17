function(kynema_setup_dependencies)
  #--------------------------------------------------------------------------
  # Required packages
  #--------------------------------------------------------------------------
  find_package(KokkosKernels REQUIRED)
  find_package(yaml-cpp REQUIRED)
  find_package(NetCDF REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(Eigen3 REQUIRED)

  #--------------------------------------------------------------------------
  # Optional packages
  #--------------------------------------------------------------------------
  #----------------------------------------
  # Sparse Direct Linear Solvers
  #----------------------------------------
  if(Kynema_ENABLE_SUPERLU)
    find_package(superlu REQUIRED)
  endif()

  if(Kynema_ENABLE_SUPERLU_MT)
    find_package(superlu_mt REQUIRED)
  endif()

  if(Kynema_ENABLE_KLU)
    find_package(KLU REQUIRED)
  endif()

  if(Kynema_ENABLE_UMFPACK)
    find_package(UMFPACK REQUIRED)
  endif()

  if(Kynema_ENABLE_MKL)
    find_package(MKL REQUIRED)
  endif()

  if(Kynema_ENABLE_CUSOLVERSP)
    if(NOT DEFINED Kokkos_ENABLE_CUDA)
      message(FATAL_ERROR "When Kynema_ENABLE_CUSOLVERSP is enabled, Kokkos must also be built with CUDA")
    endif()
  endif()

  if(Kynema_ENABLE_CUDSS)
    if(NOT DEFINED Kokkos_ENABLE_CUDA)
      message(FATAL_ERROR "When Kynema_ENABLE_CUDSS is enabled, Kokkos must also be built with CUDA")
    endif()
    find_package(cudss REQUIRED)
  endif()

  #----------------------------------------
  # GTest
  #----------------------------------------
  if(Kynema_ENABLE_TESTS)
    find_package(GTest REQUIRED)
  endif()

  #--------------------------------------------------------------------------
  # Compiler-specific settings
  #--------------------------------------------------------------------------
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
      set(FS_LIB stdc++fs)
    endif()
  endif()

  #--------------------------------------------------------------------------
  # External projects
  #--------------------------------------------------------------------------
  #----------------------------------------
  # OpenFAST/AerodynInflow (ADI) library
  #----------------------------------------
  if(Kynema_ENABLE_OPENFAST_ADI)
    find_library(OpenFast_ADI_LIBRARY NAMES aerodyn_inflow_c_binding)
    set(Kynema_ADI_LIBRARY ${OpenFast_ADI_LIBRARY} CACHE PATH "ADI library")
  endif()

  #----------------------------------------
  # ROSCO Controller library
  #----------------------------------------
  if(Kynema_ENABLE_ROSCO_CONTROLLER)
    find_library(Rosco_LIBRARY NAMES discon)
    set(Kynema_ROSCO_LIBRARY ${Rosco_LIBRARY} CACHE PATH "Rosco discon library")
  endif()
endfunction()
