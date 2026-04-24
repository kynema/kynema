include(cmake/SystemLink.cmake)

#--------------------------------------------------------------------------
# KynemaFMB Build Options
#--------------------------------------------------------------------------
macro(kynema_fmb_setup_options)
  #----------------------------------------
  # Core build options
  #----------------------------------------
  option(KYNEMA_FMB_ENABLE_TESTS "Build Tests" ON)
  option(KYNEMA_FMB_ENABLE_DOCUMENTATION "Build Documentation" OFF)
  option(KYNEMA_FMB_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  option(KYNEMA_FMB_ENABLE_IPO "Enable IPO/LTO (Interprocedural Optimization/Link Time Optimization)" OFF)
  option(KYNEMA_FMB_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

  #----------------------------------------
  # Sanitizer options
  #----------------------------------------
  option(KYNEMA_FMB_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
  option(KYNEMA_FMB_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
  option(KYNEMA_FMB_ENABLE_SANITIZER_UNDEFINED "Enable undefined behavior sanitizer" OFF)
  option(KYNEMA_FMB_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
  option(KYNEMA_FMB_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)

  #----------------------------------------
  # Build optimization options
  #----------------------------------------
  option(KYNEMA_FMB_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
  option(KYNEMA_FMB_ENABLE_PCH "Enable precompiled headers" OFF)

  #----------------------------------------
  # Static analysis options
  #----------------------------------------
  option(KYNEMA_FMB_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
  option(KYNEMA_FMB_ENABLE_CPPCHECK "Enable CppCheck analysis" OFF)

  #----------------------------------------
  # External dependencies
  #----------------------------------------
  option(KYNEMA_FMB_ENABLE_CUSOLVERSP "Use cuSolverSP for the sparse linear solver when running on CUDA Devices" OFF)
  option(KYNEMA_FMB_ENABLE_CUDSS "Use cuDSS for the sparse linear solver when running on CUDA Devices" OFF)
  option(KYNEMA_FMB_ENABLE_MKL "Use MKL for sparse linear solver when running on CPU" OFF)
  option(KYNEMA_FMB_ENABLE_KLU "Use KLU for sparse linear solver when running on CPU" OFF)
  option(KYNEMA_FMB_ENABLE_UMFPACK "Use UMFPACK for sparse linear solver when running on CPU" OFF)
  option(KYNEMA_FMB_ENABLE_SUPERLU "Use SuperLU for sparse linear solver when running on CPU" OFF)
  option(KYNEMA_FMB_ENABLE_SUPERLU_MT "Use SuperLU-MT for sparse linear solver when running on CPU" OFF)
  option(KYNEMA_FMB_ENABLE_OPENFAST_ADI "Build the OpenFAST ADI external project" OFF)
  option(KYNEMA_FMB_ENABLE_ROSCO_CONTROLLER "Build the ROSCO controller external project" OFF)
endmacro()

#--------------------------------------------------------------------------
# KynemaFMB Global Options
#--------------------------------------------------------------------------
macro(kynema_fmb_global_options)
  # Enable IPO/LTO if the option is set
  if(KYNEMA_FMB_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    kynema_fmb_enable_ipo()
  endif()
endmacro()

#--------------------------------------------------------------------------
# Project-Wide Configuration Options
#--------------------------------------------------------------------------
macro(kynema_fmb_local_options)
  #----------------------------------------
  # Core setup
  #----------------------------------------
  # Include standard project settings and create interface libraries
  include(cmake/StandardProjectSettings.cmake)
  add_library(kynema_fmb_warnings INTERFACE)
  add_library(kynema_fmb_options INTERFACE)

  #----------------------------------------
  # Compiler warnings
  #----------------------------------------
  include(cmake/CompilerWarnings.cmake)
  kynema_fmb_set_project_warnings(
    kynema_fmb_warnings
    ${KYNEMA_FMB_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
  )

  #----------------------------------------
  # Sanitizers configuration
  #----------------------------------------
  include(cmake/Sanitizers.cmake)
  kynema_fmb_enable_sanitizers(
    kynema_fmb_options
    ${KYNEMA_FMB_ENABLE_SANITIZER_ADDRESS}
    ${KYNEMA_FMB_ENABLE_SANITIZER_LEAK}
    ${KYNEMA_FMB_ENABLE_SANITIZER_UNDEFINED}
    ${KYNEMA_FMB_ENABLE_SANITIZER_THREAD}
    ${KYNEMA_FMB_ENABLE_SANITIZER_MEMORY}
  )

  #----------------------------------------
  # Build optimizations
  #----------------------------------------
  # Configure unity build
  set_target_properties(kynema_fmb_options
    PROPERTIES UNITY_BUILD ${KYNEMA_FMB_ENABLE_UNITY_BUILD}
  )

  # Configure precompiled headers
  if(KYNEMA_FMB_ENABLE_PCH)
    target_precompile_headers(
      kynema_fmb_options
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

  if(KYNEMA_FMB_ENABLE_CLANG_TIDY)
    kynema_fmb_enable_clang_tidy(kynema_fmb_options ${KYNEMA_FMB_WARNINGS_AS_ERRORS})
  endif()

  if(KYNEMA_FMB_ENABLE_CPPCHECK)
    kynema_fmb_enable_cppcheck(${KYNEMA_FMB_WARNINGS_AS_ERRORS} "")
  endif()

  #----------------------------------------
  # Coverage configuration
  #----------------------------------------
  if(KYNEMA_FMB_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    kynema_fmb_enable_coverage(kynema_fmb_options)
  endif()
endmacro()
