include(cmake/SystemLink.cmake)

#--------------------------------------------------------------------------
# Kynema Build Options
#--------------------------------------------------------------------------
macro(kynema_setup_options)
  #----------------------------------------
  # Core build options
  #----------------------------------------
  option(Kynema_ENABLE_TESTS "Build Tests" ON)
  option(Kynema_ENABLE_DOCUMENTATION "Build Documentation" OFF)
  option(Kynema_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  option(Kynema_ENABLE_IPO "Enable IPO/LTO (Interprocedural Optimization/Link Time Optimization)" OFF)
  option(Kynema_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

  #----------------------------------------
  # Sanitizer options
  #----------------------------------------
  option(Kynema_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
  option(Kynema_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
  option(Kynema_ENABLE_SANITIZER_UNDEFINED "Enable undefined behavior sanitizer" OFF)
  option(Kynema_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
  option(Kynema_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)

  #----------------------------------------
  # Build optimization options
  #----------------------------------------
  option(Kynema_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
  option(Kynema_ENABLE_PCH "Enable precompiled headers" OFF)

  #----------------------------------------
  # Static analysis options
  #----------------------------------------
  option(Kynema_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
  option(Kynema_ENABLE_CPPCHECK "Enable CppCheck analysis" OFF)

  #----------------------------------------
  # External dependencies
  #----------------------------------------
  option(Kynema_ENABLE_CUSOLVERSP "Use cuSolverSP for the sparse linear solver when running on CUDA Devices" OFF)
  option(Kynema_ENABLE_CUDSS "Use cuDSS for the sparse linear solver when running on CUDA Devices" OFF)
  option(Kynema_ENABLE_MKL "Use MKL for sparse linear solver when running on CPU" OFF)
  option(Kynema_ENABLE_KLU "Use KLU for sparse linear solver when running on CPU" OFF)
  option(Kynema_ENABLE_UMFPACK "Use UMFPACK for sparse linear solver when running on CPU" OFF)
  option(Kynema_ENABLE_SUPERLU "Use SuperLU for sparse linear solver when running on CPU" OFF)
  option(Kynema_ENABLE_SUPERLU_MT "Use SuperLU-MT for sparse linear solver when running on CPU" OFF)
  option(Kynema_ENABLE_OPENFAST_ADI "Build the OpenFAST ADI external project" OFF)
  option(Kynema_ENABLE_ROSCO_CONTROLLER "Build the ROSCO controller external project" OFF)
endmacro()

#--------------------------------------------------------------------------
# Kynema Global Options
#--------------------------------------------------------------------------
macro(kynema_global_options)
  # Enable IPO/LTO if the option is set
  if(Kynema_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    kynema_enable_ipo()
  endif()
endmacro()

#--------------------------------------------------------------------------
# Project-Wide Configuration Options
#--------------------------------------------------------------------------
macro(kynema_local_options)
  #----------------------------------------
  # Core setup
  #----------------------------------------
  # Include standard project settings and create interface libraries
  include(cmake/StandardProjectSettings.cmake)
  add_library(kynema_warnings INTERFACE)
  add_library(kynema_options INTERFACE)

  #----------------------------------------
  # Compiler warnings
  #----------------------------------------
  include(cmake/CompilerWarnings.cmake)
  kynema_set_project_warnings(
    kynema_warnings
    ${Kynema_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
  )

  #----------------------------------------
  # Sanitizers configuration
  #----------------------------------------
  include(cmake/Sanitizers.cmake)
  kynema_enable_sanitizers(
    kynema_options
    ${Kynema_ENABLE_SANITIZER_ADDRESS}
    ${Kynema_ENABLE_SANITIZER_LEAK}
    ${Kynema_ENABLE_SANITIZER_UNDEFINED}
    ${Kynema_ENABLE_SANITIZER_THREAD}
    ${Kynema_ENABLE_SANITIZER_MEMORY}
  )

  #----------------------------------------
  # Build optimizations
  #----------------------------------------
  # Configure unity build
  set_target_properties(kynema_options
    PROPERTIES UNITY_BUILD ${Kynema_ENABLE_UNITY_BUILD}
  )

  # Configure precompiled headers
  if(Kynema_ENABLE_PCH)
    target_precompile_headers(
      kynema_options
      INTERFACE
        <vector>
        <string>
        <utility>
    )
  endif()

  #----------------------------------------
  # Static analysis tools
  #----------------------------------------
  include(cmake/StaticAnalyzers.cmake)

  if(Kynema_ENABLE_CLANG_TIDY)
    kynema_enable_clang_tidy(kynema_options ${Kynema_WARNINGS_AS_ERRORS})
  endif()

  if(Kynema_ENABLE_CPPCHECK)
    kynema_enable_cppcheck(${Kynema_WARNINGS_AS_ERRORS} "")
  endif()

  #----------------------------------------
  # Coverage configuration
  #----------------------------------------
  if(Kynema_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    kynema_enable_coverage(kynema_options)
  endif()
endmacro()
