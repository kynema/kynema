function(kynema_fmb_assure_out_of_source_builds)
  # Get the real paths of the source and binary directories
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # Disallow in-source builds
  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "######################################################\n"
      "Error: In-source builds are disabled.\n"
      "Please create a separate build directory and run CMake from there.\n"
      "Example:\n"
      "  mkdir -p build\n"
      "  cd build\n"
      "  cmake ..\n"
      "######################################################"
    )
  endif()
endfunction()

# Invoke the function to enforce the out-of-source build
kynema_fmb_assure_out_of_source_builds()
