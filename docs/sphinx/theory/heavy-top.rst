.. _sec-heavy-top:

Heavy top constrained-rigid-body example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide here a simple application of the Kynema formulation for the
heavy-top problem, which is a rotating body fixed to the ground by a spherical
joint. It is a common benchmark problem for constrained-rigid-body dynamics and
for testing Lie-group time integrators like that used in Kynema. Our approach
follows much of what is described in in [@Bruls-etal:2012], but with the key
differences that we formulate the problem in inertial coordinates rather than
material coordinates and we include three translational DOFs for the spherical joint, :math:`\underline{u}_\mathrm{SJ}`.   The generalized coordinates for this problem are

.. math::

   \underline{q} = \begin{bmatrix}
    \underline{u}_\mathrm{HT} \\
    \underline{\underline{R}}_\mathrm{HT} \\
    \underline{u}_\mathrm{SJ} 
   \end{bmatrix}

where :math:`\underline{q} = \mathbb{R}^{9}`.

We assume the heavy top is a thin disk with mass
:math:`m=15` kg. The :math:`6\times6` mass matrix in material coordinates is

.. math::

   \underline{\underline{M}}^*_\mathrm{HT} = \begin{bmatrix}
   15 \mathrm{~kg}& 0 & 0 & 0 & 0 & 0\\
   0 & 15 \mathrm{~kg} & 0 & 0 & 0 & 0\\
   0 & 0 & 15 \mathrm{~kg} & 0 & 0 & 0\\
   0 & 0 & 0 & 0.234375 \mathrm{~kg~m}^2 & 0 & 0\\
   0 & 0 & 0 & 0 &  0.46875 \mathrm{~kg~m}^2 & 0\\
   0 & 0 & 0 & 0 & 0 & 0.234375 \mathrm{~kg~m}^2 \\
   \end{bmatrix}

and the :math:`9\times9` system mass matrix in material coordinates is

.. math::

   \underline{\underline{M}}^* = \begin{bmatrix}
   \underline{\underline{M}}^*_\mathrm{HT} & \underline{\underline{0}}\\
   \underline{\underline{0}} & \underline{\underline{0}}
   \end{bmatrix}

The heavy-top center of mass reference position and orientation (see
Eq. :eq:`rigid-ref`) are given by

.. math::

   \underline{x}^\mathrm{r}_\mathrm{HT} = ( 0, 1 , 0 )^T\, \mathrm{m}, \quad
   \underline{\underline{R}}^\mathrm{r}_\mathrm{HT} = \underline{\underline{I}} \,, \\

respectively. The only component of external force (see
Eq. :eq:`rbresid`) is gravity:

.. math:: \underline{f} = [0,0,-g,0,0,0,0,0,0]^T

where :math:`g=9.81` m/s\ :math:`^2`. 


The problem is constrained such that the spherical joint is located at the global origin, the heavy-top center of mass is located 1 m from the
origin, and the heavy-top body-attached-coordinate-system y-axis is pointing away from the global origin.  The Kynema implementation employs the following constraint:

.. math:: \underline{\Phi} = \begin{bmatrix}
   \underline{u}_\mathrm{HT} - \underline{u}_\mathrm{SJ} 
   - \left( \underline{\underline{R}}_\mathrm{HT} - \underline{\underline{I}}\right)
   \underline{x}^r_\mathrm{HT} \\
   \underline{u}_\mathrm{SJ}
   \end{bmatrix}

where :math:`\underline{\Phi} \in  \mathbb{R}^6`, and the contraints :math:`\underline{\lambda} \in \mathbb{R}^6` are denoted

.. math::

   \underline{\lambda} = \begin{bmatrix}
    \underline{\lambda}_\mathrm{HT} \\
    \underline{\lambda}_\mathrm{SJ} 
   \end{bmatrix}

where :math:`\underline{\lambda}_\mathrm{HT},\, \underline{\lambda}_\mathrm{SJ}
\in \mathbb{R}^3`.

The constraint gradient matrix is

.. math:: \underline{\underline{B}} = \begin{bmatrix}
     \underline{\underline{I}} & \widetilde{ \underline{\underline{R}}_\mathrm{HT} \underline{x}^{\mathrm{r}}_\mathrm{HT} } & - \underline{\underline{I}} \\
     \underline{\underline{0}} & 
     \underline{\underline{0}} & 
     \underline{\underline{I}} 
   \end{bmatrix}

and the contribution to the iteration matrix in Eq. :eq:`iteration` is

.. math:: \underline{\underline{K}}^\Phi = \begin{bmatrix}
     \underline{\underline{0}} & 
     \underline{\underline{0}} & 
     \underline{\underline{0}} \\
     \underline{\underline{0}} & 
     \widetilde{\lambda_\mathrm{HT}} 
     \widetilde{\underline{\underline{R}}_\mathrm{HT} \underline{x}^{\mathrm{r}}_\mathrm{HT} }  & \underline{\underline{0}}  \\
     \underline{\underline{0}} & 
     \underline{\underline{0}} & 
     \underline{\underline{0}} 
   \end{bmatrix}


The Kynema regression test suite includes the spinning, heavy top
problem with the following initial conditions:

.. math::

   \begin{aligned}
   \underline{u}^\mathrm{init}_\mathrm{HT} &= \left[ 0, 0, 0 \right]^T \, 
   \mathrm{m}\\
   \underline{\underline{R}}_\mathrm{HT}^\mathrm{init} &=  \underline{\underline{R}}^\mathrm{r}_\mathrm{HT} \\
   \underline{u}^\mathrm{init}_\mathrm{SJ} &= \left[ 0, 0, 0 \right]^T \, 
   \mathrm{m}
   \end{aligned}

.. math::

   \begin{aligned}
   \dot{\underline{u}}^\mathrm{init}_\mathrm{HT} &= 
   \widetilde{\omega^\mathrm{init}_\mathrm{HT}}
   \left(\underline{x}^\mathrm{r}_\mathrm{HT}
   +\underline{u}^\mathrm{init}_\mathrm{HT}\right) \\
   \omega^\mathrm{init}_\mathrm{HT} &= (0, 150, -4.61538)^T \, 
   \mathrm{rad/s}\\
   \dot{\underline{u}}^\mathrm{init}_\mathrm{SJ} &= \left[ 0, 0, 0 \right]^T\, 
   \mathrm{m/s}
   \end{aligned}

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Bruls-etal:2012

      Brüls, O., A. Cardona, and M. Arnold. 2012. “Lie Group
      Generalized-:math:`\alpha` Time Integration For Constrained
      Flexible Multibody Systems.” *Mechanism and Machine Theory*,
      121–37.
