.. _dev-plan:

Kynema-FMB Development Plan
############################

Background and overview
***********************

Kynema-FMB development started in early 2023 with primary funding from the
U.S. Department of Energy (DOE) Wind Energy Technologies Office (WETO) and with
additional funding from the DOE Exascale Computing Project (ECP). It is being
developed by researchers at the National Renewable Energy Laboratory (NREL)
and the Sandia National Laboratories (SNL).

Kynema-FMB is an open-source structural dynamics simulation code designed to
meet the research needs of WETO and the broader wind energy community for
land-based and offshore wind turbines. Kynema-FMB provides high-fidelity,
highly performant structural dynamics models that can couple with low-fidelity
aerodynamic/hydrodynamic models like those in `OpenFAST <https://github.com/OpenFAST/openfast>`_,
and high-fidelity computational fluid dynamics (CFD) models like those in the
WETO and Office of Science supported `ExaWind <https://github.com/Exawind>`_ code suite.
Kynema-FMB is designed deliberately to address shortcomings of legacy wind turbine structural
models and codes that are critical to the success of WETO modeling efforts.

Development priorities
**********************

Robustness
==========

Considering lessons learned from nearly a decade of OpenFAST development,
Kynema-FMB prioritizes software robustness through a comprehensive unit
and regression test suite, which are run through a Continuous Integration
process.  Beyond that, Kynema-FMB continuously runs a variety of static
and dynamic analysis tools to identify potential bugs.  Linters and manual
code review on every change help to ensure consistent, well designed, and
sustainable software.

Performance
===========

Kynema-FMB is performance-focused software.  Core data structures are designed
to provide optimal cache usage and all algorithms are written to best take
advantage of on-chip resources.  Kynema-FMB performs optimally on both CPU
and GPU, using hiearchical parallelism and other techniques to ensure performance
portability for problems of all sizes.

Accessibility
=============

Kynema-FMB provides a user-friendly, high level API for developers to use
to define and run their structural dynamics problems and to couple Kynema-FMB
with other codes.  This approach decouples users of Kynema-FMB from the low level
details of its implementation and improves the speed at which developers can
define and execute their problem.  Advanced users are able to use the lower level
Kynema-FMB APIs directly or to define their own interfaces, should the high
level APIs not address their needs directly.

Programming language and models
*******************************

Kynema-FMB is written in C++ with tight integration of the 
`Kokkos <https://github.com/kokkos/kokkos>`_ performance portability library.
This approach allows a single code base to achieve near-optimal performance
when run on CPU or any GPU platform, rather than requiring separate code paths.
Optimized math routines, such as those covered by BLAS and LAPACK pacakges or
sparse linear solvers, are obtained from specialized third-party libraries to ensure
state-of-the-art performance.

Key numerical algorithms
************************

Kynema-FMB models turbines using a combination of high-order nonlinear beam finite elements,
point mass elements, linear spring elements, and constraints tying them together.
For example, a turbine rotor may be modeled with three 10th-order beam elements, each
representing a blade, with their "root" nodes constrained to rotate with a hub of finite radius,
which is modeled by a point mass.  

The models necessary for mid- to high-fidelity simulation of wind turbine
structural dynamics include linear and nonlinear finite-element models coupled
through constraints equations. Kynema-FMB models These models together constitute a set of
differential-algebraic equations (DAEs) in the time domain. Kynema-FMB builds on the
experiences gained with OpenFAST, particularly its nonlinear beam-dynamics module,
`BeamDyn <https://github.com/OpenFAST/openfast/tree/main/modules/beamdyn>`_.

For more details, see the Kynema-FMB's theory documentation.

High-level development timeline
*******************************

CY = calendar year, FY = fiscal year

**CY23 Q2**: The Kynema-FMB team will implement a rigid-body dynamics solver following the
concepts described above, i.e., DAE-3 coupling, quaternion-based rotation representation, and a
generalized-alpha time integrator. This proof-of-concept implementation will be made available
in the ``main`` branch of Kynema-FMB repository and will inform the next steps in Kynema-FMB
development.

**CY23 Q3**: Implement a general GEBT-based beam element that is appropriate for constrained multi-body
simulations of a wind turbine. Enable variable order finite elements and user-defined material property
definition (appropriate for modern turbine blades). Demonstrate performance for a dynamic cantilever beam
problem and compare against `BeamDyn <https://github.com/OpenFAST/openfast/tree/main/modules/beamdyn>`_.

**CY24 Q1**: Demonstrate a wind turbine rotor simulation under prescribed loading and include code
verification results and automated testing results. Include control system
(e.g., `ROSCO <https://github.com/NREL/ROSCO/tree/main/ROSCO>`_) and pitch control of blades.
Compare simulation time against an equivalent model simulated with
`OpenFAST <https://github.com/OpenFAST/openfast>`_.

**CY24 Q3**: Demonstrate a rotor simulation with fluid-structure interaction (FSI) and a pitch control
system. Fluid will be represented in two ways. First, through a simple Blade Element Momentum Theory
(BEMT) solver and second, where the blades are represented as actuator lines in the fluid domain
(solved with the ExaWind CFD code).

**CY25 Q1**: Release a robust, well-documented, well-tested version of Kynema-FMB for land-based
turbine simulations. Demonstrate whole turbine simulation (tower, nacelle, drivetrain) capabilities
with FSI coupling to ExaWind.
