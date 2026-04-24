function(kynema_fmb_enable_coverage project_name)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Ensure that the coverage flags are set for both compilation and linking
    target_compile_options(${project_name} INTERFACE --coverage -O0 -g)
    target_link_libraries(${project_name} INTERFACE --coverage)

    # Provide a message to inform that coverage is enabled
    message(STATUS "Coverage enabled for target '${project_name}'")
  else()
    # Notify if the compiler does not support coverage
    message(STATUS "Coverage is not supported for '${CMAKE_CXX_COMPILER_ID}' compiler.")
  endif()
endfunction()
