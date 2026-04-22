# Kynema-FMB
[Documentation](https://kynema.github.io/kynema-fmb) | [Nightly test dashboard](http://my.cdash.org/index.php?project=Exawind) 

Kynema-FMB is an open-source performance portable flexible multibody (FMB) dynamics
solver designed for time-domain simulations.  While originally tailored for
wind turbine structural dynamics, the formulation and implementation are those of a general flexible-multidbody dynamics solver that can readily be applied to a wide range of systems.
Kynema was designed with a narrow focus, namely to provide a lightweight, fast,
accurate FMD solver for coupling to computational-fluid-dynamics (CFD) codes,
especially the CFD codes in the [Kynema](https://github.com/kynema) suite, for
fluid-structure-interaction (FSI) simulations.  

Kynema-FMB is equipped to model systems that can be represented as a collection of
beams and rigid bodies that are connected through constraints.  Degrees of
freedom are defined in the inertial/global frame of reference and include
displacements and rotations (formally as rotation matrices, but stored as
quaternions).  The underlying formulation is built on a Lie-group time
integrator designed for index-3 differential-algebraic equations, which is
second-order accurate in time ([Bruls et al.,
2012](https://www.sciencedirect.com/science/article/pii/S0094114X11001510)).
Beam models are based on geometrically exact beam theory and are discretized as
high-order spectral finite elements similar to those in BeamDyn ([Wang et al.,
2017](https://onlinelibrary.wiley.com/doi/full/10.1002/we.2101)).  The governing
equations for a FMD system like a wind turbine constitute a highly nonlinear
system of constrained partial-differential equations.  Kynema-FMB uses analytical
Jacobians in the nonlinear-system solves in each time step.  Linear systems use
sparse storage and several third-party sparse-linear-system solvers are
enabled. Ill conditioning of linear systems is mitigated with preconditioning
described in [Bottasso et al,
2008](https://link.springer.com/article/10.1007/s11044-007-9051-9).  Kynema-FMB is
integrated with a simple open-source controller
([ROSCO](https://github.com/NatLabRockies/ROSCO)). There is an application programming
interface (API) for coupling to geometry-resolved CFD (like that in [Sharma et
al., 2023](https://onlinelibrary.wiley.com/doi/full/10.1002/we.2886)) and
actuator-force CFD (like that in [Kuhn et al.,
2025](https://onlinelibrary.wiley.com/doi/full/10.1002/we.70010)).  In the
latter, for actuator-line models, Kynema-FMB includes an internal blade-element
solver that depends on user-provided lookup tables for coefficients of lift and
drag, i.e., aerodynamic polars.

Kynema-FMB is written in C++ and leverages Kokkos and Kokkos-Kernels
([KokkosEcosystem](https://kokkos.org/)) as its performance portability layer
enabling simulations on both CPU and GPU systems. The repository is equipped
with extensive automated testing at the unit and regression/system levels.

The following describes the high-level development objectives conceived for Kynema:
- Kynema will follow modern software development best practices, 
including test-driven development (TDD), version control,
hierarchical automated testing, and continuous integration (CI) for a
robust development environment.
- The core data structures are memory efficient and enable vectorization
and parallelization at multiple levels.
- Data structures are data-oriented to exploit methods for accelerated computing including
high utilization of chip resources (e.g., single instruction multiple data (SIMD) instruction sets) and
parallelization using GP-GPUs.
- The computational algorithms incorporate robust open-source libraries for
mathematical operations, resource allocation, and data management.
- The API design considers multiple stakeholder needs and ensure
integration with existing and future ecosystems for data science, machine learning,
and AI.
- Kynema-FMB is written in modern C++ and leverages [Kokkos](https://github.com/kokkos/kokkos)
as its performance-portability library with inspiration from the kynema stack.

## Contributing 

Kynema-FMB is an open-source project and we welcome contributions
from external developers.  To do so, open an issue describing your contribution
and make a pull request against Kynema-FMB's main branch.  Smaller contributions
are always preferred - 10 self-contained 200 line changes are easier to review
and coordinate with others than one 2000 line change.

When adding a feature, make sure that it is comprehensively covered by unit
tests and regression tests.  Bug fixes should be accompanied by at least one
test (but possibly more) which fails without the fix but now passes.

Kynema-FMB's CI process targets a number of different configurations for MacOS and
Linux, but is not fully comprehensive of the platforms we support.  In your PR,
please indicate on which platforms you've tested your contribution (i.e., Linux
x86 and CUDA v12).  This will let us know what other platforms we may have to
test against in the review process.

## Development support

Kynema-FMB is developed by researchers at The National Laboratory of the
Rockies and Sandia National Laboratories under the support of the U.S.
Department of Energy (DOE) Office of Critical Minerals and Energy Innovation
(CMEI) and the DOE Office of Science FLOWMAS Energy Earthshot Research Center.

[Kynema-FMB Software Release Record SWR-23-07](https://www.osti.gov/biblio/code-166281)

[Documentation](https://kynema.github.io/kynema-fmb/)

Send questions to michael.a.sprague@nlr.gov, Kynema Principal Investigator.
