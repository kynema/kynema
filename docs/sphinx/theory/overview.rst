.. _sec-overview:

Formulation overview
--------------------

Kynema was designed to solve flexible multibody problems that have
been discretized into a collection of :math:`k_n` nodes representing
massless points, rigid bodies, and flexible-beam nodes, where each node
has three translational degrees of freedom (DOFs) and three DOFs
defining orientation. In general, a node’s DOFs are represented on the
six-dimensional manifold, :math:`G`, with a Lie group structure,
:math:`G \in \mathbb{R}^3 \times \mathrm{SO}(3)`, where
:math:`\mathrm{SO}(3)` is the group of :math:`3\times 3` proper
orthogonal linear transformations. For example, the generalized
displacement for node :math:`i` is the pair
:math:`(\underline{u}_i,\underline{\underline{R}}_i)`, where
:math:`\underline{u}_i \in \mathbb{R}^3` is the displacement and
:math:`\underline{\underline{R}}_i\in\mathrm{SO(3)}` is the relative
rotation matrix. The associated composition operation is defined as
:math:`(\underline{u}_1,\underline{\underline{R}}_1)\circ(\underline{u}_2,\underline{\underline{R}}_2) = (\underline{u}_1+\underline{u}_2,\underline{\underline{R}}_1 \underline{\underline{R}}_2)`.
We denote the generalized displacement and velocity for node :math:`i`
as

.. math::

   [#eq:one] 

   \underline{q}_i = \begin{bmatrix}
     \underline{u}_i \\
     \underline{\underline{R}}_i
    \end{bmatrix}
   \,,\quad
    \underline{v}_i = \begin{bmatrix}
     \dot{\underline{u}}_i \\
     \underline{\omega}_i
    \end{bmatrix}
   \,,

respectively, where
:math:`\underline{q}_i \in \mathbb{R}^3\times \mathrm{SO(3)}`,
:math:`\underline{v}_i(t) \in \mathbb{R}^{6}`, an overdot denotes a time
derivative, and

.. math::

   \begin{aligned}
   \underline{\omega}_i= \mathrm{axial}\left({\dot{\underline{\underline{R}}}\, \underline{\underline{R}}^T}\right)
   \end{aligned}

is the angular velocity, for which the axial vector is defined such
that, for :math:`\underline{\underline{A}} \in \mathrm{SO(3)}` with
entries :math:`A_{ij}`,

.. math::

   \begin{aligned}
    \mathrm{axial}\left({\underline{\underline{A}}}\right) =
   \begin{bmatrix}
   A_{32}-A_{23} \\
   A_{13}-A_{31} \\
   A_{21}-A_{12}
   \end{bmatrix}
   \end{aligned}

Note that in Kynema, degrees of freedom are defined in the inertial
coordinate system.

For a discretized flexible multibody system with :math:`k_n` nodes, the
generalized displacement is organized as

.. math::

   \begin{aligned}
    \underline{q} = 
   \begin{bmatrix}
     \underline{u}_1 \\
     \underline{\underline{R}}_1 \\
     \underline{u}_2 \\
     \underline{\underline{R}}_2 \\
      \vdots \\
     \underline{u}_{k_n} \\
     \underline{\underline{R}}_{k_n} \\
    \end{bmatrix}\,, \quad
    \underline{v} = 
   \begin{bmatrix}
     \dot{\underline{u}}_1 \\
     \underline{\omega}_1 \\
     \dot{\underline{u}}_2 \\
     \underline{\omega}_2 \\
      \vdots \\
     \dot{\underline{u}}_{k_n} \\
     \underline{\omega}_{k_n} \\
    \end{bmatrix} 
   \end{aligned}

where
:math:`\underline{q} \in \left[\underline{\underline{R}}^3 \times \mathrm{SO(3)} \right]^{k_n}`,
:math:`\underline{v} \in \mathbb{R}^k`, :math:`k=6k_n`.

For a discretized flexible multibody system with :math:`k` degrees of
freedom and :math:`m` kinematic constraints, Kynema is restricted
to problems where the governing equations of motion and the constraint
equations form a residual vector
:math:`\underline{r}\in \mathbb{R}^{k+m}` where

.. math:: \underline{r} = \begin{bmatrix} \underline{R} + \underline{\underline{B}}^T \underline{\lambda} \\ \underline{\Phi} \end{bmatrix}
   :label: residual

:math:`\underline{R}\left(\underline{q},\underline{v},\dot{\underline{v}}, t \right) \in \mathbb{R}^{k}`
is the *unconstrained*-equations-of-motion residual,
:math:`\underline{\underline{B}}(\underline{q},t) \in \mathbb{R}^{m\times k}`
is the constraint-gradient matrix associated with the constraints
:math:`\underline{\Phi}(\underline{q},t)\in \mathbb{R}^m`, and
:math:`\underline{\lambda}\in \mathbb{R}^m` are the Lagrange multipliers
associated with the constraints. Kynema is restricted to problems
for which the unconstrained-equations-of-motion residual can be written
in the form

.. math:: \underline{R} = 
     \underline{\underline{M}}(\underline{q}) \dot{\underline{v}} + \underline{g}(\underline{q},\underline{v},t) - \underline{f}(t)
   %+ \uu{B}^T(\u{q},t) \u{\lambda}(t)
   :label: residual1

where :math:`\underline{\underline{M}} \in \mathbb{R}^{k\times k}` is
the mass matrix and :math:`\underline{g} \in \mathbb{R}^k` are internal
forces and :math:`\underline{f} \in \mathbb{R}^k` are external forces. The
variation of Eq. :eq:`residual` can be written

.. math:: 
   \delta \underline{r} = 
   \begin{bmatrix}
   \underline{\underline{M}}(\underline{q}) \delta \underline{\dot{v}} + \underline{\underline{G}}(\underline{q},\underline{v},t) \delta \underline{v} + \left[ \underline{\underline{K}}(\underline{q},\underline{v},\underline{\lambda},t) + \underline{\underline{K}}^\Phi(\underline{q},\underline{\lambda},t) \right] \delta \underline{q} + \underline{\underline{B}}^T \delta \underline{\lambda}\\
   \underline{\underline{B}}\, \delta \underline{q} 
   \end{bmatrix}
   :label: variation

where,
:math:`\underline{\underline{G}}, \underline{\underline{K}} \in \mathbb{R}^{k \times k}`
are the linearized damping and stiffness matrices, respectively,
:math:`\underline{\underline{K}}^\Phi \in \mathbb{R}^{k \times k}` is
the stiffness matrix associated with the constraint forces,

.. math::

   \begin{aligned}
   \delta \underline{q} = \begin{bmatrix} 
   \delta \underline{u}_1 \\
   \delta \underline{\theta}_1\\
   \delta \underline{u}_2 \\
   \delta \underline{\theta}_2\\
   \vdots \\
   \delta \underline{u}_{k_n} \\
   \delta \underline{\theta}_{k_n}\\
    \end{bmatrix}\,, \quad
   \delta \underline{v} = \begin{bmatrix} 
   \delta \dot{\underline{u}}_1 \\
   \delta \underline{\omega}_1\\
   \delta \dot{\underline{u}}_2 \\
   \delta \underline{\omega}_2\\
   \vdots \\
   \delta \dot{\underline{u}}_{k_n} \\
   \delta \underline{\omega}_{k_n}\\
    \end{bmatrix}
   \end{aligned}

:math:`\delta \underline{q}, \delta \underline{v} \in \mathbb{R}^k`, and
:math:`\delta \underline{u}_i, \delta \underline{\theta}_i \in \mathbb{R}^3`
are the virtual displacement and virtual rotation, respectively, in
inertial coordinates associated with node :math:`i`.

In the following, we describe the time-integration algorithm for index-3
differential-algebraic-equation (DAE-3) systems, which is the backbone
of the Kynema framework. We then discuss the governing equations
for a single rigid body, and then the theory and numerical
discretization for flexible beams and practical consideration in
modeling wind turbine blades. The full constrained system for a
land-based turbine is described, including the interface for
fluid-structure-interaction simulations where the fluid is simulated by
an external solver. We finish the formulation with a description of the
computational implementation.
