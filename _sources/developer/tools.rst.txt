Developer Tools
===============

This page describes static and dynamic analysis tools used in the development
of Kynema-FMB. These tools are run and must pass as part of the CI process,
so it will likely become important for integration into a developer's process.

clang-format
------------

`ClangFormat <https://clang.llvm.org/docs/ClangFormat.html>`_ is used for
linting to enforce a consistent code style. It can be installed with most package
managers.

ClangFormat is configured by the ``.clang-format`` file at the top of the repository.
If the tool is run from the top directory, it will automatically detect and load the
settings in the configuration file.

clang-tidy
----------

`ClangTidy <https://clang.llvm.org/extra/clang-tidy/>`_ is another linting tool
tool which enforces a variety of rules on the code in order to avoid common
bugs.

ClangTidy is configure by the ``.clang-tidy`` file at the top of the repository.
To run it, configure Kynema-FMB with the ``KYNEMA_FMB_ENABLE_CLANG_TIDY`` option.

Cppcheck
--------

`Cppcehck <https://cppcheck.sourceforge.io/>`_ is yet another linting tools
which detects undefined behavior and dangerous constructs with very few
false positives.

To run Cppcheck, configure Kynema-FMB with the ``KYNEMA_FMB_ENABLE_CPPCHECK`` option.
