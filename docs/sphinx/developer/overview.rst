.. _overview:

Overview
========

Kynema-FMB is developed on GitHub usng git for source control and the GoogleTest
framework for unit and regression testing.  Becoming fluent with these tools
will make the development experience much smoother.  What follows is an overview
of Kynema-FMB's structure and expectations for contributers.

.. contents::

Version control
---------------

Kynema-FMB follows a process of continuous deployment with a single ``main`` branch
which for which all tests must always pass.  Therefore, in order for changes to be
incorperated into the main branch, they should be tested and reviewed thoroughly.
This requirement extends beyond simply passing all CI checks, but it should also
include adding new tests to cover newly added code and running on all availible platforms.

Kynema-FMB periodically creates stable, versioned release branches.  These branches
will only ever receive bug fixes and should be used when API and solution stability
are more important than access to the latest features or performance improvements.

Pull requests
-------------

When adding a new feature, open the pull request (PR) on GitHub against the ``main``
branch.  These PRs should contain a description of the code to be added as well as
a link to any GitHub issues that are addressed.  A list of platforms where the code
has been tested and a discussion of any performance changes or anticipated limitations
should also be included.

PRs will require a code review from a core developer to ensure that the code meets quality,
robutness, and performance standards.  Reviewers should be prepared to provide prompt and
detailed reviews regarding both code design and its fitness for inclusion within Kynema-FMB.

PRs should be kept small, both in terms of lines of code and conceptually.  A series of
ten 200 line PRs implementing a change will likely be easier to read and review than a
single 2000 line PR.  Similarly, just performing one change, even if it has to proliferate
to a hundred different places will help with the review process.  These rules will also
help to avoid conflicts with other developers and to cultivate a useful git history.

GitHub issues
-------------

Prior to starting work on a new feature or fixing a bug, it is helpful to start a new
issue on GitHub.  This issue will help to coordinate with other developers as well
as to discuss design decisions early in the development process.  A PR which has
been previously documented and discussed in an issue will be much easier to review
and more likely to be accepted than one submitted without any preface.

Testing
-------

All pull requests should be maximally tested in the test suite. This
includes new code as well as existing code that is modified. In particular, all bug
fixes must be accompanied by tests which are broken on the ``main`` branch but
are fixed by the PR.  

Kynema-FMB includes both unit tests and regression tests.  While the distinction
between these is not always obvious, in general a unit test will test a single
behavior of a single part of the code in a controlled environment with an expected,
analytical result.  In comparison, a regression test is typically more "end-to-end"
in nature, involving many parts of the code, with the the test focusing on ensuring
that the values produce have not been changed by code changes rather than on absolute
correctness.  Both types of test are important and exercise the code differently, so
new features should ideally be accompanied by both, when it is reasonable to do so.

Documentation
-------------

Kynema-FMB has three sources of documentation: this guide, Doxygen documentation, and
documentation tests.

When changing the behavior of Kynema-FMB, it's important to keep this document up to date.
The theory section in particular should reflect the implementation, as this will provide
confidence and reporducibility to users and other developers.  Other changes, such as
chaging compilation options or code restructuring should also be reflected here to ensure
consistency.

Functions and classes added to Kynema-FMB should be properly annotated to provide Doxygen
documentation.  When editing existing code, it is easy for the Doxygen documentation to fall
out of date, so be pay careful attention to keeping everything in sync.

Kynema-FMB also provides documentation tests - these are stand-alone programs which act as tutorials,
with extensive comments for users explaining the logic.  When adding a new user-facing feature,
add a documentation test explaining how to set up a problem using it and what facets a user
may need to consider when adapting it to their own needs.  Note that these documentation tests
are also important for testing our installation process, so they will need to be added to that
CI script, which will also ensure that said test does not fall out of date.
