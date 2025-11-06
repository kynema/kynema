.. _sec-heavy-top:

Heavy top constrained-rigid-body example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide here a simple application of the Kynema formulation for
the heavy-top problem, which is a rotating body fixed to the ground
by a spherical joint. It is a common benchmark problem for
constrained-rigid-body dynamics and for testing Lie-group time
integrators like that used in Kynema. We follow the problem
description found in [@Bruls-etal:2012], but with the key difference
that we formulate the problem in inertial coordinates rather than
material coordinates. We assume the heavy top is a thin disk with mass
:math:`m=15` kg. The :math:`6\times6` mass matrix in material
coordinates is

.. math::

   \underline{\underline{M}}^* = \begin{bmatrix}
   15 \mathrm{~kg}& 0 & 0 & 0 & 0 & 0\\
   0 & 15 \mathrm{~kg} & 0 & 0 & 0 & 0\\
   0 & 0 & 15 \mathrm{~kg} & 0 & 0 & 0\\
   0 & 0 & 0 & 0.234375 \mathrm{~kg~m}^2 & 0 & 0\\
   0 & 0 & 0 & 0 &  0.46875 \mathrm{~kg~m}^2 & 0\\
   0 & 0 & 0 & 0 & 0 & 0.234375 \mathrm{~kg~m}^2 \\
   \end{bmatrix}

The heavy-top center of mass reference position and orientation (see
Eq. :eq:`rigid-ref`) are given by

.. math::

   \underline{x}^\mathrm{r} = ( 0, 1 , 0 )^T\, \mathrm{m}, \quad
   \underline{\underline{R}}^\mathrm{r} = \underline{\underline{I}} \,, \\

respectively. The only component of external force (see
Eq. :eq:`rbresid`) is gravity:

.. math:: \underline{f} = [0,0,-g,0,0,0]^T

where :math:`g=9.81` m/s\ :math:`^2`. The problem is constrained such
that the center of mass is located 1 m from the origin, which can be
written as three constraint equations as

.. math:: \underline{\Phi} = \underline{\underline{R}}\, \underline{x}^\mathrm{r} - \underline{x}^c \in  \mathbb{R}^3

where :math:`\underline{\Phi} \in \mathbb{R}^3`, :math:`\underline{x}^c` is the current center-of-mass position,
and for which the constraint gradient matrix is

.. math:: \underline{\underline{B}}  = \begin{bmatrix}
   -\underline{\underline{I}} & \widetilde{- \underline{\underline{R}}\, \underline{x}^\mathrm{r}}
   \end{bmatrix}

:math:`\underline{\underline{B}} \in \mathbb{R}^{3 \times 6}`. The
stiffness matrix associated with linearization of the constraint forces (see Eq. :eq:`variation`) is

.. math:: \underline{\underline{K}}^\Phi = \begin{bmatrix} 
   \underline{\underline{0}} & \underline{\underline{0}}\\
   \underline{\underline{0}} & 
   \widetilde{\lambda} \, \widetilde{\underline{\underline{R}} \underline{x}^\mathrm{r}}
   \end{bmatrix}

where :math:`\underline{\lambda} \in  \mathbb{R}^3` are the Lagrange multipliers.  The Kynema regression test suite includes the spinning, heavy top
problem with the following initial conditions:

.. math::

   \begin{aligned}
   \underline{u}^\mathrm{init} &= \left[ 0, 0, 0 \right]^T \, \mathrm{m}\\
   \underline{\underline{R}}^\mathrm{init} &=  \underline{\underline{R}}^\mathrm{r}
   \end{aligned}

.. math::

   \begin{aligned}
   \omega^\mathrm{init} &= (0, 150, -4.61538)^T \, \mathrm{rad/s}\\
   \dot{\underline{u}}^\mathrm{init} &= \widetilde{\omega^\mathrm{init}}\left(\underline{x}^\mathrm{r}+\underline{u}^\mathrm{init}\right)\, \mathrm{m/s}
   \end{aligned}

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Bruls-etal:2012

      Brüls, O., A. Cardona, and M. Arnold. 2012. “Lie Group
      Generalized-:math:`\alpha` Time Integration For Constrained
      Flexible Multibody Systems.” *Mechanism and Machine Theory*,
      121–37.
