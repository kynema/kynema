function(kynema_set_project_warnings project_name WARNINGS_AS_ERRORS CLANG_WARNINGS GCC_WARNINGS CUDA_WARNINGS)

    # Set default Clang warnings if not provided
    if("${CLANG_WARNINGS}" STREQUAL "")
        set(CLANG_WARNINGS
      -Wall               # Enable all standard warnings
      -Wextra             # Enable extra warnings
      -Wshadow            # Warn if a variable declaration shadows one from a parent context
      -Wnon-virtual-dtor  # Warn if a class with virtual functions has a non-virtual destructor
      -Wcast-align        # Warn for potential performance problem casts
      -Wunused            # Warn on anything being unused
      -Woverloaded-virtual # Warn if you overload (not override) a virtual function
      -Wpedantic          # Warn if non-standard C++ is used
      -Wconversion        # Warn on type conversions that may lose data
      -Wsign-conversion   # Warn on sign conversions
      # -Wnull-dereference  # Warn if a null dereference is detected (disabled due to false positive in g++)
      -Wdouble-promotion  # Warn if float is implicitly promoted to double
      -Wformat=2          # Warn on security issues around functions that format output (e.g., printf)
      -Wimplicit-fallthrough # Warn on statements that fall through without an explicit annotation
    )
    endif()

    # Set default GCC warnings if not provided
    if("${GCC_WARNINGS}" STREQUAL "")
        set(GCC_WARNINGS
      ${CLANG_WARNINGS}         # Include Clang warnings
      -Wmisleading-indentation  # Warn if indentation implies blocks where blocks do not exist
      -Wduplicated-cond         # Warn if an if/else chain has duplicated conditions
      -Wduplicated-branches     # Warn if if/else branches have duplicated code
      -Wlogical-op              # Warn about logical operations being used where bitwise were probably wanted
      -Wuseless-cast            # Warn if you perform a cast to the same type
      -Wsuggest-override        # Warn if an overridden member function is not marked 'override' or 'final'
    )
    endif()

    # Set default CUDA warnings if not provided
    if("${CUDA_WARNINGS}" STREQUAL "")
        set(CUDA_WARNINGS
      -Wall
      -Wextra
      -Wunused
      -Wconversion
      -Wshadow
    )
    endif()

    # Treat warnings as errors if specified
    if(WARNINGS_AS_ERRORS)
        message(STATUS "Warnings are treated as errors for project '${project_name}'")
        list(APPEND CLANG_WARNINGS -Werror)
        list(APPEND GCC_WARNINGS -Werror)
        list(APPEND CUDA_WARNINGS -Werror)
    endif()

    # Set compiler-specific warning flags
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(PROJECT_WARNINGS_CXX ${CLANG_WARNINGS})
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(PROJECT_WARNINGS_CXX ${GCC_WARNINGS})
    else()
        message(WARNING "No compiler warnings set for CXX compiler: '${CMAKE_CXX_COMPILER_ID}'")
    endif()

    set(PROJECT_WARNINGS_CUDA "${CUDA_WARNINGS}")

    # Apply the warning flags to the target
    target_compile_options(
    ${project_name}
    INTERFACE
    $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_WARNINGS_CXX}>
    $<$<COMPILE_LANGUAGE:CUDA>:${PROJECT_WARNINGS_CUDA}>
  )

endfunction()
