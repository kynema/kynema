.. _sec-rb-springs:

Rigid body with three springs
-----------------------------

This section provides model details for a rigid-body solver with three
massless taut-line "mooring lines", represented by geometrically
nonlinear springs. The model configuration is appropriate for modeling
the DeepCWind floating platform.

The generalized model degrees of freedom are

.. math::

   \underline{q} = \begin{bmatrix} 
     \underline{u} \\
     \underline{\underline{R}} \\
     \underline{u}^\mathrm{sp1} \\
     \underline{u}^\mathrm{sp2} \\
     \underline{u}^\mathrm{sp3}
    \end{bmatrix} \quad
    \underline{v} = \begin{bmatrix} 
     \dot{\underline{u}} \\
     \underline{\omega} \\
     \dot{\underline{u}}^\mathrm{sp1} \\
     \dot{\underline{u}}^\mathrm{sp2} \\
     \dot{\underline{u}}^\mathrm{sp3}
    \end{bmatrix}

:math:`\underline{u} \in \mathbb{R}^3` and
:math:`\underline{\underline{R}} \in \mathbb{R}^{3\times 3}` are rigid
body displacement and orientation, respectively,
:math:`\underline{u}^{\mathrm{sp}i} \in \mathbb{R}^6` is spring i end
displacement (2 nodes with 3 DOF each; see :ref:`sec-spring`) and

.. math::

   \underline{u}^{\mathrm{sp}i} = 
   \begin{bmatrix}
   \underline{u}^{\mathrm{sp}i}_1 \\ \underline{u}^{\mathrm{sp}i}_2 
   \end{bmatrix}

The system input parameters include the following:

- :math:`\underline{\underline{M}}^* \in \mathbb{R}^{6\times 6}`: Rigid
  body mass matrix in material coordinate system

- :math:`\underline{x}^\mathrm{r} \in \mathbb{R}^3`: Rigid body reference
  position

- :math:`\underline{\underline{R}}^\mathrm{r} \in \mathbb{R}^{3\times3}`: Rigid
  body reference orientation

- :math:`\underline{u}^\mathrm{init} \in \mathbb{R}^3`: Rigid body
  initial displacement

- :math:`\underline{\underline{R}}^\mathrm{init} \in \mathbb{R}^{3\times3}`:
  Rigid body initial orientation

- :math:`\dot{\underline{u}}^\mathrm{init} \in \mathbb{R}^3`: Rigid body
  initial velocity

- :math:`\underline{\omega}^\mathrm{init} \in \mathbb{R}^3`: Rigid body
  initial angular velocity

- :math:`\underline{x}^{\mathrm{sp,r}i} \in \mathbb{R}^6`: Spring
  :math:`i` reference location

- :math:`k^{\mathrm{sp}i} \in \mathbb{R}`: Spring :math:`i` spring
  constant

- :math:`L^{\mathrm{sp,r}i} \in \mathbb{R}`: Spring :math:`i` unstretched
  length

and where

.. math::

   \underline{x}^{\mathrm{sp,r}i} =
   \begin{bmatrix}
   \underline{x}^{\mathrm{sp,r}i}_1 \\ \underline{x}^{\mathrm{sp,r}i}_2
   \end{bmatrix}

The constraints for each spring can be written

.. math::

   \begin{aligned}
    \underline{u}_1^{\mathrm{sp}i} &=  \underline{u} + \underline{\underline{R}} \underline{r}^{\mathrm{sp,r}i}_1
   - \underline{r}^{\mathrm{sp,r}i}  = \underline{u} + \left(\underline{\underline{R}}-\underline{\underline{I}}\right) 
   \underline{r}^{\mathrm{sp,r}i}_1 \\
   \underline{u}_2^{\mathrm{sp}i} &= 0
   \end{aligned}

where :math:`\underline{r}^{\mathrm{sp,r}i}_1 
= \underline{x}^{\mathrm{sp,r}i}_1 - \underline{x}^\mathrm{r}`.

The full-system residual required for time integration can be written as follows (see Eq. :eq:`residual`):

.. math::

   \underline{R} = \begin{bmatrix}
   \underline{\underline{M}}\, \dot{\underline{v}} +\underline{g} - \underline{f} + \underline{\underline{B}}^T\underline{\lambda} \\
   \underline{\Phi}
   \end{bmatrix} 
   \in \mathbb{R}^{42}

where :math:`\underline{\lambda} \in \mathbb{R}^{18}`, :math:`\underline{\underline{B}} \in \mathbb{R}^{18 \times 24}`

.. math::

   \underline{\underline{M}} 
   = \begin{bmatrix}
   \underline{\underline{M}}^\mathrm{rb} & \underline{\underline{0}}_{6\times 18}\\
   \underline{\underline{0}}_{18\times 6} & \underline{\underline{0}}_{18\times 18}
   \end{bmatrix} \in \mathbb{R}^{24\times 24}
   \,,\quad \underline{g} = \begin{bmatrix} 
   \underline{g}^\mathrm{rb} \\ 
   \underline{g}^\mathrm{sp1} \\
   \underline{g}^\mathrm{sp2} \\
   \underline{g}^\mathrm{sp3} 
   \end{bmatrix} \in \mathbb{R}^{24}
   \,,\quad
   \underline{f}
   = 
   \begin{bmatrix}
   \underline{f}^\mathrm{rb}\\
   \underline{0}_{18}\\
   \end{bmatrix} 
   \in \mathbb{R}^{24}

In the above, :math:`\underline{\underline{M}}^\mathrm{rb} \in \mathbb{R}^{6 \times 6}` and :math:`\underline{g}^\mathrm{rb}\in \mathbb{R}^6` are defined in :ref:`sec-rigid`, :math:`\underline{f}^\mathrm{rb}\in \mathbb{R}^6` is the force applied to the rigid-body center of mass, and :math:`\underline{g}^{\mathrm{sp}i}\in\mathbb{R}^6` is defined in :ref:`sec-spring`. The constraints are given by

.. math::

   \underline{\Phi} =
   \begin{bmatrix}
   \underline{u}_1^{\mathrm{sp}1} - \underline{u} 
   - \left(\underline{\underline{R}}-\underline{\underline{I}} \right) \underline{r}^{\mathrm{sp,r}1} \\
   \underline{u}_2^{\mathrm{sp}1} \\
   \underline{u}_1^{\mathrm{sp}2} - \underline{u} 
   - \left(\underline{\underline{R}}-\underline{\underline{I}} \right) \underline{r}^{\mathrm{sp,r}2} \\
   \underline{u}_2^{\mathrm{sp}2} \\
   \underline{u}_1^{\mathrm{sp}3} - \underline{u} 
   - \left(\underline{\underline{R}}-\underline{\underline{I}} \right) \underline{r}^{\mathrm{sp,r}3} \\
   \underline{u}_2^{\mathrm{sp}3} \\
   \end{bmatrix} \quad
   \in \mathbb{R}^{18}

The specific entries of :math:`\underline{\underline{B}}` are given by

.. math::

   \begin{aligned}
   \underline{\underline{B}} =
   \begin{bmatrix}
   -\underline{\underline{I}}_{3 \times 3} & \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r}1} } 
   & \underline{\underline{I}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} 
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} \\
   %
   \underline{\underline{0}}_{3 \times 3} &  \underline{\underline{0}}_{3 \times 3}
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{I}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} 
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} \\
   %
   -\underline{\underline{I}}_{3 \times 3} & \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r}2} } 
   & \underline{\underline{0}}_{3 \times 3}  & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{I}}_{3 \times 3} 
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} \\
   %
   \underline{\underline{0}}_{3 \times 3} &  \underline{\underline{0}}_{3 \times 3}
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} 
   & \underline{\underline{I}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} \\
   %
   -\underline{\underline{I}}_{3 \times 3} & \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r}3} } 
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3}
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{I}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3}\\
   %
   \underline{\underline{0}}_{3 \times 3} &  \underline{\underline{0}}_{3 \times 3}
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} 
   & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{0}}_{3 \times 3} & \underline{\underline{I}}_{3 \times 3} \\
   \end{bmatrix}
   \end{aligned}

The full matrices :math:`\underline{\underline{C}}`, :math:`\underline{\underline{K}}`, and :math:`\underline{\underline{K}}^\Phi` for the iteration matrix Eq. :eq:`iteration` in :ref:`sec-gen-alpha` are given by

.. math::

   \underline{\underline{C}} = \begin{bmatrix}
   \underline{\underline{C}}^\mathrm{rb} & \underline{\underline{0}}_{6 \times 18}\\
   \underline{\underline{0}}_{18 \times 6} & \underline{\underline{0}}_{18 \times 18}
   \end{bmatrix} 
   \in \mathbb{R}^{24 \times 24}

.. math::

   \underline{\underline{K}} = \begin{bmatrix}
   \underline{\underline{K}}^\mathrm{rb} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6}\\
   \underline{\underline{0}}_{6 \times 6} & \underline{\underline{K}}^\mathrm{sp1} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6}\\
   \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{K}}^\mathrm{sp2} & \underline{\underline{0}}_{6 \times 12}\\
   \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6} &  \underline{\underline{0}}_{6 \times 6} & \underline{\underline{K}}^\mathrm{sp3}
   \end{bmatrix}
   \in \mathbb{R}^{24 \times 24}

.. math::

   \underline{\underline{K}}^\Phi = \begin{bmatrix}
   \underline{\underline{0}}_{3 \times 3} & 
   \underline{\underline{0}}_{3 \times 3} & 
   \underline{\underline{0}}_{3 \times 18} \\
   \underline{\underline{0}}_{3 \times 3} & 
   \left( 
    \widetilde{\lambda}_2  \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r1}}} +  
    \widetilde{\lambda}_4  \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r2}}} +  
    \widetilde{\lambda}_6  \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r3}}}\right)  &
   \underline{\underline{0}}_{3 \times 18} \\
   \underline{\underline{0}}_{18 \times 3} & 
   \underline{\underline{0}}_{18 \times 3} & 
   \underline{\underline{0}}_{18 \times 18} 
   \end{bmatrix} \in \mathbb{R}^{24 \times 24}

.. note::

   The term :math:`\underline{\underline{K}}^\Phi` is include here for
   completeness, but is not in the current Kynema implementation, and
   nonlinear-system convergence has been satisfactory.


