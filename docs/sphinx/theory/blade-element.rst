.. _sec-blade-element:

Blade-element aerodynamics solver
---------------------------------

In this section we describe the theory for a simple blade-element solver
for fluid-structure interaction calculations.

For beams coupled to an aerodynamic solver, we follow the WindIO
requirement that the aerodynamic reference line is the same as the
structure reference line.

.. note::
 
   This section is under construction and needs to be finalized to match
   the Kynema implementation.


Initialization
~~~~~~~~~~~~~~

Our blade-element solver requires the following user inputs for each
beam:

- :math:`\eta^\mathrm{ac}_j \in [0,1]`,
  :math:`j\in\{1,\ldots n^\mathrm{ac}\}`, where
  :math:`\eta^\mathrm{ac}_j` are the nondimensional locations of the
  aerodynamic sections along the reference line and
  :math:`n^\mathrm{ac}` is the number of aerodynamic sections.

- :math:`\tau^\mathrm{ac}_j` is the section aerodynamic twist

- :math:`c^\mathrm{ac}_j` is the section chord length

- location of the aerodynamic center in the local aerodynamic section
  plane with respect to the reference line; :math:`y_j^\mathrm{ac}`,
  :math:`z_j^\mathrm{ac}`

- :math:`C_j = C_j(\alpha_j)` which is a function that gives
  :math:`C^L_j`, :math:`C^D_j`, :math:`C^M_j`, the coefficients of lift,
  drag, and moment per unit span, respecitvely.

Calculate initialization quantities required by the FSI API (see
§\ `[sec:fsi] <#sec:fsi>`__)

.. math::

   \begin{aligned}
   n^\mathrm{motion} = n^\mathrm{force} = n^\mathrm{ac} \\
   \xi_j^\mathrm{motion,map}=\xi_j^\mathrm{force,map} =2 \eta^\mathrm{ac}_j-1
   \end{aligned}

:math:`\underline{x}_j^\mathrm{motion,map,0}` and
:math:`\widehat{q}_j^\mathrm{motion,map,0}` are calculated in the FSI
API based on :math:`\xi_j^\mathrm{motion,map}` and the underlying basis
functions.

.. math::

   \begin{aligned}
   \underline{x}_j^\mathrm{motion,0} = \underline{x}_j^\mathrm{force,0} = \underline{x}_j^\mathrm{motion,map,0}
   + \underline{\underline{R}}\left( \widehat{q}^\mathrm{motion,map,0}_j\right)  
   \begin{bmatrix} 0 \\ 
   y_j^\mathrm{ac} \\
   z_j^\mathrm{ac} 
   \end{bmatrix}
   \end{aligned}

Calculate quantities required by the blade-element solver

.. math::

   \begin{aligned}
     \Delta s_j =  \left \{
   \begin{array}{ll}
       \int_{\xi_1^\mathrm{motion,map}}^{(\xi_j^\mathrm{motion,map}+\xi_{j+1}^\mathrm{motion,map})/2} J(\xi) d \xi 
       &  j =  1 \\
       \int_{(\xi_{j-1}+\xi_j)/2}^{(\xi_{j}+\xi_{j+1})/2} J(\xi) d \xi
       &  j =  \{2,3,\ldots,n^\mathrm{motion}-1\} \\
       \int_{(\xi_{j-1}+\xi_{j})/2}^{\xi^\mathrm{motion,map}_{j}} J(\xi) d \xi 
       &  j =  n^\mathrm{motion}
   \end{array}
   \right .
   \end{aligned}

Force calculations based on blade-element polars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each aerodynamic section, we must calculate the relative velocity at
the aerodynamic center based on the inflow velocity at that point and
the velocity of the aerodynamic center. That relative velocity is
projected onto the aerodynamic section plane and is used to calculate
lift, drag, and moment forces.

(1) For a given inflow velocity
:math:`\dot{\underline{u}}_j^\mathrm{inflow}` and velocity of the
aerodynamic center :math:`\dot{\underline{u}}_j^\mathrm{fl}`, calculate
the relative flow velocity at point
:math:`\underline{x}_j^\mathrm{motion}`:

.. math::

   \begin{aligned}
     \dot{\underline{u}}^\mathrm{rel} =
   \begin{bmatrix}
   0 & 0 & 0 \\
   0 & 1 & 0 \\
   0 & 0 & 1
   \end{bmatrix}
   \left[ \underline{\underline{R}}\left(\hat{q}_j^\mathrm{motion,map}\right) \underline{\underline{R}}\left( \widehat{q}_j^\mathrm{motion,map,0}\right) \right]^T 
   \left( \dot{\underline{u}}^\mathrm{inflow}_j- \dot{\underline{u}}^\mathrm{fl}_j \right)
   \end{aligned}

(2) Given aerodynamic twist :math:`\tau_j`, calculate the angle of
attack as :math:`\alpha_j = \beta_j - \tau_j`, where

.. math::

   \begin{equation}
   \beta_j = \left\{
   \begin{array}{ll}
   \mathrm{arccos} \left( \frac{\dot{\underline{u}}^\mathrm{rel}_j \cdot \hat{i}_y}{| \dot{\underline{u}}^\mathrm{rel}_j |}\right) & \mathrm{if}\, \dot{\underline{u}}^\mathrm{rel}_j \cdot \hat{i}_y \ge 0 \\
   2 \pi - \mathrm{arccos} \left( \frac{\dot{\underline{u}}^\mathrm{rel}_j \cdot \hat{i}_y}{| \dot{\underline{u}}^\mathrm{rel}_j |}\right) & \mathrm{if}\, \dot{\underline{u}}^\mathrm{rel}_j \cdot \hat{i}_y < 0 
   \end{array}
   \right .
   \end{equation}

(3) Calculate :math:`C^L_j`, :math:`C^D_j`, :math:`C^M_j` given
:math:`\alpha_j` and calculate the force and moment in the aerodynamic
coordinates:

.. math::

   \begin{aligned}
     \underline{f}_j = 
   \begin{bmatrix}
   0 \\
   \left( C^D_j \cos \tau_j - C^L_j \sin \tau_j   \right) \\
   \left( C^D_j \sin \tau_j + C^L_j \cos \tau_j   \right) 
   \end{bmatrix}
   \frac{1}{2} \rho c_j \Delta s_j |\dot{\underline{u}}^\mathrm{rel}_j|^2 
   \end{aligned}

.. math::

   \begin{aligned}
     \underline{m}_j = 
   \begin{bmatrix}
   C^M_j  \\
   0 \\
   0 
   \end{bmatrix}
   \frac{1}{2} \rho c^2_j \Delta s_j |\dot{\underline{u}}^\mathrm{rel}_j|^2 
   \end{aligned}

(4) Calculate force and moment in intertial coordiates (see
Eqs. (`[eq:force] <#eq:force>`__)-(`[eq:moment] <#eq:moment>`__)):

.. math::

   \begin{aligned}
   \underline{f}^\mathrm{force}_j = \underline{\underline{R}}\left(\widehat{q}_j^\mathrm{motion,map} \right) 
   \underline{\underline{R}}\left(\widehat{q}_j^\mathrm{motion,map,0}\right) \underline{f}_j\\
   \underline{m}^\mathrm{force}_j = \underline{\underline{R}}\left(\widehat{q}_j^\mathrm{motion,map} \right) 
   \underline{\underline{R}}\left(\widehat{q}_j^\mathrm{motion,map,0}\right) \underline{m}_j
   \end{aligned}

