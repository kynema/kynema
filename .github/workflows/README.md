# GitHub Workflow Definitions for CI

This collection of GitHub workflow definitions is designed for Continuous Integration (CI) and is organized as follows:

- **correctness-linux.yaml**: Builds and runs the entire test suite on `ubuntu-latest` using both _gcc_ and _clang_ compilers. Both `Debug` and `Release` builds are tested, with `AddressSanitizer`, `LeakSanitizer`, and `UndefinedBehaviorSanitizer` enabled.  The CppCheck and ClangTidy static analysis tools are both enabled for all configurations.

- **correctness-macos.yaml**: Builds and runs the entire test suite on `macos-latest` using the _clang_ compiler, both the default version provided by Apple and the latest version. Similar to the Linux workflow, both `Debug` and `Release` builds are tested, with `AddressSanitizer` and `UndefinedBehaviorSanitizer` enabled.  The CppCheck and ClangTidy static analysis tools are both enabled for all configurations.

- **deploy-docs.yaml**: Builds the Sphinx and Doxygen documentation.  When run on the main branch, this also deploys the documentation so that it is accessible online.

- **formatting.yaml**: Checks formatting for all source, header, and test files. Note that this workflow does not automatically fix formatting errors; developers are responsible for correcting any issues detected.

- **install-spack.yaml**: Builds and installs Kynema-FMB using the Spack package manager and then builds each of the documentation tests against this installation.  This script not only tests the documentation tests themselves, also ensures that our Spack package and CMake installation scripts are working properly for our end users.
