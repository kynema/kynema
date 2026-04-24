# Function to enable various sanitizers for a given target
function(kynema_fmb_enable_sanitizers
  project_name
  ENABLE_SANITIZER_ADDRESS
  ENABLE_SANITIZER_LEAK
  ENABLE_SANITIZER_UNDEFINED_BEHAVIOR
  ENABLE_SANITIZER_THREAD
  ENABLE_SANITIZER_MEMORY)

  # Initialize an empty list to hold the enabled sanitizers
  set(SANITIZERS "")

  # Check if using GCC or Clang and enable the appropriate sanitizers
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")

    if(ENABLE_SANITIZER_ADDRESS)
      list(APPEND SANITIZERS "address")
    endif()

    if(ENABLE_SANITIZER_LEAK)
      list(APPEND SANITIZERS "leak")
    endif()

    if(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR)
      list(APPEND SANITIZERS "undefined")
    endif()

    if(ENABLE_SANITIZER_THREAD)
      if("address" IN_LIST SANITIZERS OR "leak" IN_LIST SANITIZERS)
        message(WARNING "Thread sanitizer cannot be used with Address or Leak sanitizers.")
      else()
        list(APPEND SANITIZERS "thread")
      endif()
    endif()

    if(ENABLE_SANITIZER_MEMORY AND CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
      message(WARNING "Memory sanitizer requires MSan-instrumented code, otherwise it may report false positives.")
      if("address" IN_LIST SANITIZERS OR "thread" IN_LIST SANITIZERS OR "leak" IN_LIST SANITIZERS)
        message(WARNING "Memory sanitizer cannot be used with Address, Thread, or Leak sanitizers.")
      else()
        list(APPEND SANITIZERS "memory")
      endif()
    endif()

  # Check if using MSVC and enable the appropriate sanitizers
  elseif(MSVC)
    if(ENABLE_SANITIZER_ADDRESS)
      list(APPEND SANITIZERS "address")
    endif()

    if(ENABLE_SANITIZER_LEAK OR ENABLE_SANITIZER_UNDEFINED_BEHAVIOR OR ENABLE_SANITIZER_THREAD OR ENABLE_SANITIZER_MEMORY)
      message(WARNING "MSVC only supports the Address sanitizer.")
    endif()
  endif()

  # Join the sanitizers into a comma-separated string
  list(JOIN SANITIZERS "," LIST_OF_SANITIZERS)

  # If any sanitizers are enabled, apply them to the target
  if(LIST_OF_SANITIZERS)
    if(NOT MSVC)
      # Apply sanitizers for GCC/Clang.
      target_compile_options(${project_name} INTERFACE -fsanitize=${LIST_OF_SANITIZERS})
      target_link_options(${project_name} INTERFACE -fsanitize=${LIST_OF_SANITIZERS})
    else()
      # Apply sanitizers for MSVC.
      string(FIND "$ENV{PATH}" "$ENV{VSINSTALLDIR}" index_of_vs_install_dir)
      if("${index_of_vs_install_dir}" STREQUAL "-1")
        message(
          SEND_ERROR
          "Using MSVC sanitizers requires setting the MSVC environment before building the project. Please manually open the MSVC command prompt and rebuild the project."
        )
      endif()
      target_compile_options(${project_name} INTERFACE /fsanitize=${LIST_OF_SANITIZERS} /Zi /INCREMENTAL:NO)
      target_compile_definitions(${project_name} INTERFACE _DISABLE_VECTOR_ANNOTATION _DISABLE_STRING_ANNOTATION)
      target_link_options(${project_name} INTERFACE /INCREMENTAL:NO)
    endif()
  endif()

endfunction()
