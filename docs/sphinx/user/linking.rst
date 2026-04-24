Linking To Kynema-FMB
======================

Kynema-FMB is fully discoverable using CMake's ``find_package`` utility.
Simply add the following line to your ``CMakeLists.txt`` file.

.. code-block:: cmake
   
    find_package(Kynema-FMB REQUIRED)

This utility will search your path for your Kynema-FMB installation and load its target information.
To link against Kynema-FMB, add to your ``CMakeLists.txt`` the line

.. code-block:: cmake

   target_link_libraries(my_executable PRIVATE Kynema-FMBFMB::kynema_fmb_library)

This line will link to Kynema-FMB and all of its dependencies - there is no need to include any transitive dependencies, such as Kokkos, explicitly.
