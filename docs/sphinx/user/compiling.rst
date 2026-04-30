Compiling
=========

Kynema-FMB is developed in C++20 and is designed to be buildable on any system with a compliant compiler.
It utilizes `Kokkos <https://github.com/kokkos/kokkos>`_ to ensure performance portability, allowing it to run on any platform supported by these projects.
We strive to test Kynema-FMB on a wide range of platforms, including Linux and macOS, although it is not feasible to cover every possible configuration.
This document outlines the build procedure verified to work on Linux (RHEL8).
For additional assistance tailored to your specific setup, please contact the developers.

Spack Installation
------------------

The easiest way to use Kynema-FMB is through the `Spack <https://spack.io/>`_ package manager.
Once you have downloaded and set up Spack for your environment, simply run

.. code-block:: bash

    spack install kynema-fmb

To see the latest list of supported configuration options, check out the package file or run

.. code-block:: bash

    spack info kynema-fmb

Once it is installed, you can load the Kynema-FMB library and its dependencies into your environment using

.. code-block:: bash

    spack load kynema-fmb

Development using Spack Developer Workflow
------------------------------------------

One easy way to set up a development environment for Kynema-FMB is to use Spack's Developer Workflow.
To setup an environment for working on Kynema-FMB, setup Spack and then run the following commands:

.. code-block:: bash

    mkdir kynema-fmb
    cd kynema-fmb
    spack env create -d .
    spack env activate .
    spack add kynema-fmb+tests
    spack install
    spack develop kynema-fmb@main
    spack concretize -f
    spack install

Kynema-FMB's source code will now be located in the kynema-fmb folder, but can be accessed from anywhere by

.. code-block:: bash

    spack cd -c kynema-fmb

After editing the code here, it can be rebuilt by running

.. code-block:: bash

    spack install

To run the tests, first access the build folder through the spack command

.. code-block:: bash

    spack cd -b kynema-fmb

Next, the tests can be run either through ctest or directly from the unit test or regression test executables

.. code-block:: bash

    ctest
    ./tests/unit_tests/kynema_fmb_unit_tests
    ./tests/regression_tests/kynema_fmb_regression_tests

You can also build Kynema-FMB from this folder using standard make commands.

For more information, please see Spack's documentation:
https://spack-tutorial.readthedocs.io/en/latest/tutorial_developer_workflows.html

Building and Developing in Kynema-FMB Directly
-----------------------------------------------

The following sections outline how to build and develop Kynema-FMB without Spack's Developer Workflows.
The main complication here is that developers will have to manage their environment and dependencies manually, which may be an unnecessary complication or a freeing feature, depending on your perspective.

Dependencies
------------

Before building Kynema-FMB, you'll need the following:

- C++ compiler that supports the C++20 standard
- `CMake <https://cmake.org/>`_: the default build system for C++ projects, version 3.21 or later
- `Kokkos <https://github.com/kokkos/kokkos>`_: core programming model for performance portability
- `KokkosKernels <https://github.com/kokkos/kokkoskernels>`_: performance portable linear algebra library
- `netCDF <https://github.com/Unidata/netcdf-c>`_: I/O data Format
- `Suite-Sparse <https://github.com/DrTimothyAldenDavis/SuiteSparse>`_: For the KLU sparse direct solver.  Other solvers, such as SuperLU are also possible to use.
- `Eigen <http://eigen.tuxfamily.org/>`_: For the Eigen linear algebra package
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_: A package for reading YAML files, to be used by regression tests
- `GoogleTest <https://github.com/google/googletest>`_: unit testing package

Installing Third Party Libraries
--------------------------------

There are several methods to obtain the necessary Third Party Libraries (TPLs) for building Kynema-FMB, however the simplest is to use the `spack <https://spack.io/>`_ package manager.
Spack offers a comprehensive set of features for development and dependency management.
The following is a quick-start guide for installing and loading the TPLs required to build Kynema-FMB.

Clone the spack repository, load the spack environment, and let spack learn about your system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    git clone git@github.com:spack/spack.git
    source spack/share/spack/setup-env.sh
    spack compiler find
    spack external find

Install GoogleTest, netCDF, Suite-Sparse, and Eigen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    spack install googletest
    spack install netcdf-c
    spack install eigen
    spack install suite-sparse

Install Kokkos and Kokkos Kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a simple serial build

.. code-block:: bash

    spack install kokkos
    spack install kokkos-kernels


To compile with OpenMP support for parallelism on CPU based machines

.. code-block:: bash

    spack install kokkos+openmp
    spack install kokkos-kernels+openmp

To compile with CUDA support

.. code-block:: bash

    spack install kokkos+cuda+wrapper
    spack install kokkos-kernels+cuda+cublas

To compile with ROCm support

.. code-block:: bash

    spack install kokkos+rocm
    spack install kokkos-kernels+rocblas

Load the TPLs into your environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    spack load googletest
    spack load suite-sparse
    spack load netcdf-c
    spack load eigen
    spack load kokkos
    spack load kokkos-kernels

Building Kynema-FMB
--------------------

The following is written assuming the TPLs in hand and the environment configured as described above.

Clone Kynema-FMB and setup a build directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    git clone git@github.com:kynema/kynema-fmb.git
    cd kynema-fmb
    mkdir build
    cd build

Configure cmake
~~~~~~~~~~~~~~~

When building Kynema-FMB, you must specify which sparse direct solver package you want to use.
We support many options here, but the default recommendation is to use suite-sparse's KLU solver for CPU builds.

For a CPU-based build which includes building unit tests, you can configure with KLU using the command

.. code-block:: bash

    cmake ../ -DKynema-FMB_ENABLE_KLU=ON

If Kokkos was built with CUDA support, you will need to use the nvcc_wrapper for compilation.
You will also get your choice of native CUDA solvers (CUDSS or cuSolverSP).
For best performance, CUDSS is currently recommended.

.. code-block:: bash

    cmake ../ -DCMAKE_CXX_COMPILER=nvcc_wrapper -DKynema-FMB_ENABLE_CUDSS=ON

You can also use any CPU-based direct solver with a CUDA build.
You may want to do this to reduce memory usage on device, or it may be faster for your problem.
In this case, the system matrix and residual are calculated on GPU, copied to host for the solve step, and then the solution is copied back to GPU.
For this mode of operation, simply configure Kynema-FMB as

.. code-block:: bash

    cmake ../ -DCMAKE_CXX_COMPILER=nvcc_wrapper -DKynema-FMB_ENABLE_KLU=ON

If Kokkos was built with ROCm support, you will need to use the hipcc program for compilation.
Currently, we do not support any native solvers for ROCm, so a CPU based solver (such as KLU) must be used.

.. code-block:: bash

    cmake ../ -DCMAKE_CXX_COMPILER=hipcc -DKynema-FMB_ENABLE_KLU=ON

Build and Test
~~~~~~~~~~~~~~

Currently, Kynema-FMB builds several shared libraries by default.
To ensure that their unit tests pass, these libraries must be copied into the directory where the tests are executed.

.. code-block:: bash

    make -j
    ctest --output-on-failure

Once built, the unit test executable can also be run directly from the build directory

.. code-block:: bash

    ./tests/unit_tests/kynema_fmb_unit_tests

External Controllers
~~~~~~~~~~~~~~~~~~~~

Wind turbine simulations often use shared library controllers to prescribe movements to blades.
While Kynema-FMB supports calling to any shared library provided by the user, it can also detect the ROSCO controller if it is in the system path.
To turn on this feature, configure CMake with the command

.. code-block:: bash

    cmake ../ -DKynema-FMB_ENABLE_ROSCO_CONTROLLER=ON

This option will define the convenience global variable `Kynema-FMB_ROSCO_LIBRARY`, which is a string containing the location of the ROSCO library and can be used to initiaize Kynema-FMB's controller wrapper.

Similarly, Kynema-FMB can call to OpenFAST's AeroDyn module as a shared library to provide an aerodynamic inflow model.
To find this library, if it is in the system path, configure Kynema-FMB with the command

.. code-block:: bash

    cmake ../ -DKynema-FMB_ENABLE_OPENFAST_ADI=ON

This option will define the convenience global variable `Kynema-FMB_ADI_LIBRARY`, which is a string containing the location of the AeroDyn library, which can be used to initialize Kynema-FMB's AeroDyn inflow wrapper.

Build Options
-------------

Kynema-FMB has several build options which can be set either when running
CMake from the command line or through a GUI such as ccmake.

- ``Kynema-FMB_ENABLE_CLANG_TIDY`` enables the Clang-Tidy static analysis tool
- ``Kynema-FMB_ENABLE_COVERAGE`` enables code coverage analysis using gcov
- ``Kynema-FMB_ENABLE_CPPCHECK`` enables the CppCheck static analysis tool
- ``Kynema-FMB_ENABLE_IPO`` enables link time optimization
- ``Kynema-FMB_ENABLE_PCH`` builds precompiled headers to potentially decrease compilation time
- ``Kynema-FMB_ENABLE_SANITIZER_ADDRESS`` enables the address sanitizer runtime analysis tool
- ``Kynema-FMB_ENABLE_SANITIZER_LEAK`` enables the leak sanitizer runtime analysis tool
- ``Kynema-FMB_ENABLE_SANITIZER_MEMORY`` enables the memory sanitizer runtime analysis tool
- ``Kynema-FMB_ENABLE_SANITIZER_THREAD`` enables the thread sanitizer runtime analysis tool
- ``Kynema-FMB_ENABLE_SANITIZER_UNDEFINED`` enables the undefined behavior sanitizer runtime analysis tool
- ``Kynema-FMB_ENABLE_TESTS`` builds Kynema-FMB's test suite
- ``Kynema-FMB_ENABLE_UNITY_BUILD`` uses unity builds to potentially decrease compilation time
- ``Kynema-FMB_WRITE_OUTPUTS`` builds Kynema-FMB with VTK support for visualization in tests. Will need the VTK TPL to be properly configured
- ``Kynema-FMB_WARNINGS_AS_ERRORS`` treats warnings as errors, including warnings from static analysis tools
- ``Kynema-FMB_ENABLE_KLU`` builds Kynema-FMB with support for Suite-Sparse's KLU solver; in our experience, this is solver is fast and robust for many of our problems.
- ``Kynema-FMB_ENABLE_UMFPACK`` builds Kynema-FMB with support for Suite-Sparse's UMFPACK solver.
- ``Kynema-FMB_ENABLE_SUPERLU`` builds Kynema-FMB with support forthe  SuperLU solver
- ``Kynema-FMB_ENABLE_SUPERLU_MT`` builds Kynema-FMB with support for SuperLU-mt, a threaded version of SuperLU which may be configured to run in parallel on CPU.
- ``Kynema-FMB_ENABLE_MKL`` builds Kynema-FMB with MKL's sparse direct solver, which can take advantage
  of multiple threads to run in parallel on CPU.
- ``Kynema-FMB_ENABLE_CUDSS`` builds Kynema-FMB with cuDSS, the next generation sparse direct solver of CUDA; still in pre-release at the time of writing, it is the preferred CUDA based solver if the platform supports it.
- ``Kynema-FMB_ENABLE_CUSOLVERSP`` builds Kynema-FMB with the cuSolverSP sparse direct solver.
- ``Kynema-FMB_ENABLE_ROSCO_CONTROLLER`` detects the ROSCO controller shared library and defines the `Kynema-FMB_ROSCO_LIBRARY` variable
- ``Kynema-FMB_ENABLE_OPENFAST_ADI`` detects the OpenFAST AeroDyn shared library and defines the `Kynema-FMB_ADI_LIBRARY` variable
