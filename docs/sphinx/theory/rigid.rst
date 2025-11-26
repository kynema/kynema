.. _sec-rigid:

Rigid body
----------

In this section we describe the terms required to calculate the residual
vector and iteration matrix in Algorithm 1 for simulation of the
dynamics of a rigid body with six degrees of freedom. A Kynema
rigid body has reference position and orientation given as

.. math:: \underline{q}^\mathrm{r} = \begin{bmatrix} 
     \underline{x}^\mathrm{r} \\
     \underline{\underline{R}}^\mathrm{r} \\
    \end{bmatrix} 
    :label: rigid-ref

where :math:`\underline{q}^\mathrm{r} \in \mathbb{R}^3\times \mathrm{SO(3)}`. The
generalized degrees of freedom (displacement and rotation) are given by

.. math::

   \begin{aligned}
    \underline{q} = \begin{bmatrix} 
     \underline{u} \\
     \underline{\underline{R}} \\
    \end{bmatrix} 
   \end{aligned}

where :math:`\underline{q} \in \mathbb{R}^3\times \mathrm{SO(3)}`, such
that current position and orientation, in inertial coordinates, are
given by

.. math::

   \begin{aligned}
    \underline{x}^\mathrm{c} = \underline{x}^\mathrm{r} + \underline{u}\\
    \underline{\underline{R}}^\mathrm{c} = \underline{\underline{R}}\,\underline{\underline{R}}^\mathrm{r}
   \end{aligned}

respectively. The rigid body is defined by a mass matrix defined in
material coordinates, and is notated as
:math:`\underline{\underline{M}}^* \in \mathbb{R}^{6\times 6}`, where
the asterisk superscript denotes material-coordiate definition

The unconstrained governing equations can be written in residual form,
following Eq. :eq:`residual1`, as

.. math:: \underline{R} = \underline{\underline{M}}\, \dot{\underline{v}} +\underline{g} - \underline{f} 
   :label: rbresid

where
:math:`\underline{R}, \underline{g}, \underline{f} \in\mathbb{R}^6`,

.. math::
   \underline{\underline{M}} 
   = \begin{bmatrix}
   m \underline{\underline{I}} & m \tilde{\eta}^T \\
   m \tilde{\eta} & \underline{\underline{\rho}}
   \end{bmatrix} 
   = \underline{\underline{\mathcal{RR}^\mathrm{r}}}\, \underline{\underline{M}}^* {\underline{\underline{\mathcal{RR}^\mathrm{r}}}}^T \in \mathbb{R}^{6\times6}

is the mass matrix in intertial coordinates, which defines :math:`m`,
:math:`\underline{\eta}\in \mathbb{R}^3`,
:math:`\underline{\underline{\rho}}\in\mathbb{R}^{3\times 3}`,

.. math::

   \underline{g} = \begin{bmatrix}
   m \tilde{\omega} \tilde{\omega} \underline{\eta}  \\ 
   \tilde{\omega} \underline{\underline{\rho}} \underline{\omega} 
   \end{bmatrix}

and where

.. math::

   \underline{\underline{\mathcal{RR}^\mathrm{r}}}=
   \begin{bmatrix}
   \underline{\underline{R}}~\underline{\underline{R}}^\mathrm{r}& \underline{\underline{\mathrm{r}}} \\
   \underline{\underline{R}} & \underline{\underline{R}}~\underline{\underline{R}}^\mathrm{r}
   \end{bmatrix}
   \in \mathbb{R}^{6\times6}

Variation of Eq. :eq:`rbresid` can be written as

.. math::

   \begin{aligned}
   \delta \underline{R} = \underline{\underline{M}} \delta \underline{\dot{v}} + \underline{\underline{G}} \delta \underline{v} + \underline{\underline{K}} \delta \underline{q} 
   \label{eq:rbvariation}
   \end{aligned}

where

.. math::

   \begin{aligned}
   \underline{\underline{G}} =
   \begin{bmatrix}
   \underline{\underline{r}} & \widetilde{ \widetilde{\omega} m \underline{\eta} }^T
            + \widetilde{\omega} m \widetilde{\eta}^T\\
   \underline{\underline{\mathrm{r}}} & \widetilde{\omega} \underline{\underline{\rho}} - \widetilde{\underline{\underline{\rho}} \underline{\omega}}
   \end{bmatrix}
   \end{aligned}

.. math::

   \begin{aligned}
   \underline{\underline{K}} =
   \begin{bmatrix}
   \underline{\underline{0}} & \left( \dot{\widetilde{\omega}} + \tilde{\omega}\tilde{\omega}
           \right) m \widetilde{\eta}^T\\
   \underline{\underline{0}} & \ddot{\widetilde{u}} m \widetilde{\eta}
            + \left(\underline{\underline{\rho}}\dot{\widetilde{\omega}}
                    -\widetilde{\underline{\underline{\rho}} \dot{\underline{\omega}}} \right)
            + \widetilde{\omega} \left( \underline{\underline{\rho}} \widetilde{\omega}
            - \widetilde{ \underline{\underline{\rho}}\underline{\omega}} \right)
   \end{bmatrix}
   \end{aligned}

For a multibody system, the rigid-body residual is inserted into the
appropriate rows of the global residual,
Eq. :eq:`residual`, and matrices
:math:`\underline{\underline{M}}`, :math:`\underline{\underline{G}}`,
and :math:`\underline{\underline{K}}` are assembled, via direct
stiffness summation, into their global counterparts in the iteration
matrix, Eq. :eq:`iteration`.


