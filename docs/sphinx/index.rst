===========
Kynema
===========

Overview
========

`Kynema <https://github.com/kynema/kynema>`_ [@Sprague-Mudafort:2022,@Sprague-etal:2025]
is an open-source flexible multibody dynamics (FMD) solver designed for time-domain simulations.  While originally created for wind turbine structural dynamics, the formulation and implementation are those of a general FMD solver and can readily be applied to other systems.  Kynema was designed with a narrow focus, namely to provide a lightweight, accurate FMD solver for coupling to fluid-dynamics codes in wind turbine research, especially the `ExaWind <https://github.com/exawind>`_ [@Sprague-etal:2020,@Sharma-etal:2023,@Kuhn-etal:2025] suite of computational-fluid-dynamics codes.
Wind turbine blades and towers are long slender structures; as such turbines can be represented at high-fidelity with beams, rigid bodies, and constraints.  Kynema provides these model elements, where degrees of freedom are defined in the inertial/global frame of reference and include displacements and rotations (formally as rotation matrices, but stored as quaternions).
The underlying formulation is built on a Lie-group time integrator for index-3 differential-algebraic equations which is second-order accurate in time [@Bruls-etal:2012].
Beam models are based on geometrically exact beam theory (GEBT) and are discretized as high-order spectral finite elements similar to those in the BeamDyn module [@Wang-etal:2017] of `OpenFAST <https://github.com/openfast/>`_ [@Jonkman:2013].
The governing equations for a FMD system like a wind turbine form a highly nonlinear system of constrained partial-differential equations.  Kynema uses analytical Jacobians in the nonlinear-system solves performed at each time step.  Linear systems use sparse storage and several third-party sparse-linear-system solvers are enabled. Ill conditioning of the linear systems are mitigated with preconditioning described in [@Bottasso-etal:2008].  Kynema is integrated with a simple open-source controller [@Abbas-etal:2022].
Kynema is written in C++ and leverages Kokkos and Kokkos-Kernels as its performance portability layer enabling simulations on both CPU and GPU systems.   The repository is equipped with extensive automated testing at the unit and regression/system levels including several full reference megawatt-scale reference turbines.
Kynema fills the need for a lightweight, open-source turbine structural dynamics code that is high fidelity, robust, fast, and capable on running on different computer architectures.

Software-development objectives of Kynema
==============================================

- Kynema adheres to modern software development best practices. The
  development process emphasizes test-driven development (TDD), version control,
  hierarchical automated testing, and continuous integration, leading to a robust
  development environment.
- Kynema is developed in modern C++ and leverages
  `Kokkos <https://github.com/kokkos/kokkos>`_ as its performance-portability
  library, drawing inspiration from the ExaWind stack.
- The core data structures are crafted to be memory efficient, enabling
  vectorization and parallelization at multiple levels.
- These structures are data-oriented to leverage accelerated computing methods,
  including high utilization of chip resources (e.g., single instruction multiple
  data, SIMD), parallelization through GPGPUs or other hardware, and support for
  memory-efficient architectures.
- The computational algorithms incorporate robust open-source libraries for
  mathematical operations, resource allocation, and data management.
- The API design considers the needs of multiple stakeholders, ensuring
  integration with existing and future ecosystems for data science, machine
  learning, and AI.

Table of contents
=================

.. toctree::
   :maxdepth: 2

   user/user
   theory/index
   developer/index
   testing/index
   publications
   acknowledgement

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

References (this page)
======================

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Abbas-etal:2022

      Abbas, N.J., D.S. Zalkind, L. Pao, and A. Wright.
      2022. 
      “A reference open-source controller for fixed and floating 
      offshore wind turbines."
      *Wind Energy Science* **7**
      53-73.
      https://doi.org/10.5194/wes-7-53-2022

   .. container:: csl-entry
      :name: ref-Bauchau:2011

      Bauchau, O. A. 2011. *Flexible Multibody Dynamics*. Springer.

   .. container:: csl-entry
      :name: ref-Bottasso-etal:2008

      Bottasso, C.L., D. Dopicao, and L. Trainelli. 2008. “On the 
      optimal scaling of index three {DAEs} in multibody dynamics."
      *Multibody System Dynamics* **19**
      3--20.
      https://doi.org/10.1007/s11044-007-9051-9

   .. container:: csl-entry
      :name: ref-Bruls-etal:2012

      Brüls, O., A. Cardona, and M. Arnold. 2012. “Lie Group
      Generalized-:math:`\alpha` time integration for constrained
      flexible multibody systems.” *Mechanism and Machine Theory* **48**,
      121–37.
      https://doi.org/10.1016/j.mechmachtheory.2011.07.017

   .. container:: csl-entry
      :name: ref-Jonkman:2013

      Jonkman, J. M. 2013. “The new modularization framework for the
      FAST wind turbine CAE tool.” In *Proceedings of the 51st AIAA
      Aerospace Sciences Meeting Including the New Horizons Forum and
      Aerospace Exposition*. Grapevine, Texas.
      https://www.osti.gov/servlets/purl/1068607 

   .. container:: csl-entry
      :name: ref-Kuhn-etal:2025

      Kuhn, M., M. Henry de Frahan, and P. Mohan et al. 2025. “AMR-Wind:
      A Performance-Portable, High-Fidelity Flow Solver for Wind Farm
      Simulations.” *Wind Energy* **28**, e70010.
      https://doi.org/https://doi.org/10.1002/we.70010

   .. container:: csl-entry
      :name: ref-Sharma-etal:2023

      Sharma, A., M. J. Brazell, and G. Vijayakumar et al. 2023.
      “ExaWind: Open-Source CFD for Hybrid-RANS/LES Geometry-Resolved
      Wind Turbine Simulations in Atmospheric Flows.” *Wind Energy* **27**
      (3): 225–57. https://doi.org/10.1002/we.2886

   .. container:: csl-entry
      :name: ref-Sprague-etal:2020

      Sprague, M.A., S. Ananthan, G. Vijayakumar, and M. Robinson. 2020.
      "ExaWind: A multi-fidelity modeling and simulation environment for 
      wind energy." *Journal of Physics: Conference Series* **1452**, 012071.
      https://doi.org/10.1088/1742-6596/1452/1/012071

   .. container:: csl-entry
      :name: ref-Sprague-Mudafort:2022

      Sprague, M.A. and R. Mudafort, 2022, Kynema [SWR-23-07].
      https://doi.org/10.11578/dc.20251006.3

   .. container:: csl-entry
      :name: ref-Sprague-etal:2025

      Sprague, M., D. Slaughter, F. Bhuiyan, D. Dement, and P. Crozier, 2025,
      Introducing Kynema, an open-source performance-portable
      flexible-multibody-dynamics solver. Presented at NAWEA/WindTech 2025,
      Richardson, TX, October 15-17.  Published as Technical report
      NREL/PR-5000-97597.  https://docs.nrel.gov/docs/fy26osti/97597.pdf

   .. container:: csl-entry
      :name: ref-Wang-etal:2017

      Wang, Q., M. A. Sprague, J. Jonkman, N. Johnson, and B. Jonkman.
      2017. “BeamDyn: A High-Fidelity Wind Turbine Blade Solver in the
      FAST Modular Framework.” *Wind Energy* **20**, 1439–62.
      https://doi.org/10.1002/we.2101 

