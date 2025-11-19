.. _sec-lsfe:

Legendre spectral finite elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The GEBT equations are discretized with Legendre spectral finite
elements (LSFEs), which are defined by :math:`P` nodes and
:math:`(P-1)`\ th order Lagrangian-interpolant basis functions. Nodes
are located at the :math:`P` Gauss-Legendre-Lobatto (GLL) points (in the
element reference domain), :math:`\xi_i` for
:math:`i\in\{1,2, \ldots,P\}`, which are the solutions to the equation

.. math::

   \begin{aligned}
    \left( 1-\xi^2\right) \frac{\partial L_{P-1}}{\partial \xi} = 0
   \end{aligned}

for :math:`\xi \in[-1,1]`, where :math:`L_{P-1}(\xi)` is the Legendre
polynomial of degree :math:`(P-1)` [@Deville-etal:2002]. The beam
reference line definition is determined by the nodal values

.. math::

   \begin{aligned}
    \underline{q}^\mathrm{r}_i = \begin{bmatrix}
     \underline{x}^\mathrm{r}_i \\
     \hat{q}^\mathrm{r}_i \\
   \end{bmatrix}
   \end{aligned}

where :math:`\underline{x}^\mathrm{r}_i \in \mathbb{R}^3` is position and
:math:`\hat{q}^\mathrm{r}_i \in \mathbb{R}^4` is the quanternion representation of the reference orientation matrix :math:`\underline{\underline{R}}^\mathrm{r}_i \in \mathrm{SO}(3)`.


The beam reference line position is represented as 

.. math::

   \begin{aligned}
   \underline{x}^{\mathrm{r},h}(\xi) =  \sum_{\ell=1}^{P} \phi_\ell(\xi) \underline{x}^\mathrm{r}_\ell
   \end{aligned}

and the quaternion orientation representation is 

.. math::

   \begin{aligned}
   \widehat{q}^\mathrm{r} (\xi) =
   \frac{\sum_{\ell=1}^P \widehat{q}^\mathrm{r}_\ell
    \phi_\ell(\xi) }
   {\left \Vert \sum_{\ell=1}^P \widehat{q}^\mathrm{r}_\ell
    \phi_\ell (\xi) \right \Vert }
   \end{aligned}


where :math:`\phi_\ell(\xi)` is the Lagrangian-interpolant of the
:math:`\ell^{th}` node. For LSFEs, those can be written as

.. math::

   \begin{aligned}
   \phi_\ell(\xi) = \frac{-1}{P (P-1)} 
   \frac{\left(1-\xi^2\right) \frac{\partial L_{P-1}}{\partial \xi}}{(\xi - \xi_\ell) L_{P-1}(\xi_\ell)}
   \end{aligned}

for :math:`\xi \in [-1,1]`. 

The orientation of the beam reference line at any :math:`\xi` is

.. math::

   \begin{aligned}
   \underline{\underline{R}}^r\left[\hat{q}^r\left(\xi\right)\right] 
   \end{aligned}

where, for any quaternion 

.. math::
   \widehat{q} = \begin{bmatrix} q \\ \underline{q} \end{bmatrix}

(:math:`q \in \mathbb{R}`, :math:`\underline{q} \in \mathbb{R}^3`, :math:`\left \Vert \widehat{q} \right \Vert = 1`)

.. math::

   \begin{aligned}
   \underline{\underline{R}}\left(\hat{q}\right) = 
	    \underline{\underline{I}} + 2 q \widetilde{q} + 2 \widetilde{q} \widetilde{q}
	   \end{aligned}


The weak form of the residual of the
governing equations Eq. :eq:`stronggoverning`,
for the :math:`i^\mathrm{th}`-node, can be written (after integration by
parts) as

.. math::
   \underline{R}_i = 
   \int_{-1}^{1} \phi_i \underline{\mathcal{R}}\, J d\xi = 
   \int_{-1}^{1} \left[
   \frac{\partial \phi_i}{\partial \xi} 
   \left(\underline{\mathcal{F}}^\mathrm{E1} +
   \underline{\mathcal{F}}^\mathrm{D1} \right)J^{-1} 
   + \phi_i \left(
   \underline{\mathcal{F}}^\mathrm{I}  
   + \underline{\mathcal{F}}^\mathrm{E2}  
   + \underline{\mathcal{F}}^\mathrm{D2}  
   - \underline{\mathcal{F}}^\mathrm{ext}
   \right)
   \right]
    J d\xi
   :label: weakresidual

:math:`\underline{R}_i\in\mathbb{R}^6`, for all
:math:`i\in\{1,\ldots,P\}`, where we have mapped the domain
:math:`s\in[0,L]` to :math:`\xi \in[-1,1]`, and
:math:`J(\xi) \in \mathbb{R}` is the Jacobian of the mapping,

.. math::

   \begin{aligned}
   J(\xi) = \sqrt{\frac{\partial \underline{x}^{\mathrm{r},h}}{\partial \xi}^T \frac{\partial \underline{x}^{\mathrm{r},h}}{\partial \xi} }
   \end{aligned}

We remark that standard linear (i.e., 2-node) and quadratic (i.e.,
3-node) elements are a subset of LSFEs.

As described above, the generalized displacement, :math:`\underline{q}`
is in :math:`\mathbb{R}^3 \times \mathrm{SO(3)}`. However, rotations are
represented as quaternions for storage and interpolation, and are stored
at nodes. In this form, nodal degrees of freedom are denoted, for the
:math:`i^\mathrm{th}` node,

.. math::

   \begin{aligned}
    \underline{q}_i = \begin{bmatrix}
     \underline{u}_i \\
     \hat{q}_i \\
   \end{bmatrix}
   \end{aligned}

where :math:`\underline{u}_i \in \mathbb{R}^3` and
:math:`\hat{q}_i \in \mathbb{R}^4`. The generalized displacement in
:math:`\mathbb{R}^3\times\mathrm{SO(3)}` along the beam reference line
is then given as

.. math::

   \begin{aligned}
    \underline{q}(s,t) = \begin{bmatrix} \underline{u}^h \\ \underline{\underline{R}} \left( \widehat{q}^h \right)
   \end{bmatrix}
   \end{aligned}

where displacement is interpolated in the normal manner, i.e.,

.. math::

   \begin{aligned}
   \underline{u}^h =  \sum_{j=1}^{p} \phi_j\underline{u}_j
   \end{aligned}

but quaternion interpolation requires normalization, i.e.,

.. math::

   \begin{aligned}
   \widehat{q}^h = \frac{ \sum_{j=1}^{p} \phi_j \hat{q}_j }
   {|| \sum_{j=1}^{p} \phi_j \hat{q}_j ||}
   \end{aligned}


Introducing a quadrature scheme with :math:`n^Q` points with locations
and weights, :math:`\xi_k^Q`, :math:`w_k^Q`,
:math:`k\in \{1, \ldots, n^Q\}`, respectively, the approximate form of
the residual, Eq. :eq:`weakresidual`, can be
written

.. math::

   \begin{aligned}
   \underline{R} = \begin{bmatrix}
   \underline{F}^\mathrm{I}_1 + \underline{F}^\mathrm{E}_1 
   + \underline{F}^\mathrm{D}_1 
   - \underline{F}^\mathrm{ext}_1  \\
   \underline{F}^\mathrm{I}_2 + \underline{F}^\mathrm{E}_2 
   + \underline{F}^\mathrm{D}_2 
   - \underline{F}^\mathrm{ext}_2  \\
   \vdots \\
   \underline{F}^\mathrm{I}_P + \underline{F}^\mathrm{E}_P 
   + \underline{F}^\mathrm{D}_P 
   - \underline{F}^\mathrm{ext}_P 
   \end{bmatrix}
   \end{aligned}

where :math:`\underline{R} \in \mathbb{R}^{6 P}` and

.. math::

   \begin{aligned}
   \underline{F}_i^\mathrm{I} &=
   \sum_{k=1}^{n^Q}
   J(\xi^Q_k) \phi_i(\xi^Q_k) \underline{\mathcal{F}}^I(\xi^Q_k) w^Q_k\, \\
   \underline{F}_i^\mathrm{E} &=
   \sum_{k=1}^{n^Q}
   \left[\phi_i^\prime(\xi^Q_k) 
   {\underline{\mathcal{F}}^{E1}}(\xi^Q_k)+ J(\xi^Q_k) \phi_i(\xi^Q_k) 
   \underline{\mathcal{F}}^\mathrm{E2}(\xi^Q_k) \right] w^Q_k\, \\
   \underline{F}_i^\mathrm{D} &=
   \sum_{k=1}^{n^Q}
   \left[\phi_i^\prime(\xi^Q_k) 
   {\underline{\mathcal{F}}^{D1}}(\xi^Q_k)+ J(\xi^Q_k) \phi_i(\xi^Q_k) 
   \underline{\mathcal{F}}^\mathrm{D2}(\xi^Q_k) \right] w^Q_k\, \\
   \underline{F}_i^\mathrm{ext} &=
   \sum_{k=1}^{n^Q} \phi_i (\xi^Q_k)
   \underline{F}^{ext}(\xi^Q_k) J(\xi^Q_k) w^Q_k 
   \end{aligned}

and 

.. math::
    \phi_i^\prime(\xi^Q_k) = \left .\frac{\partial \phi_i}{\partial \xi}\right|_{\xi=\xi^Q_k}


The matrices required for the time-integration iteration matrix in
Eq. :eq:`iteration` are constructed from the
following:

.. math::

   \begin{aligned}
   \underline{\underline{M}}_{ij} =
   \sum_{k=1}^{n^Q} &
   \phi_i(\xi^Q_k) \underline{\underline{M}}(\xi^Q_k) \phi_j(\xi^Q_k) J(\xi^Q_k) w^Q_k \\
   %
   \underline{\underline{G}}_{ij} =
   \sum_{k=1}^{n^Q} &
   \Big\{
   \phi'_i(\xi^Q_k) \underline{\underline{D}}(\xi^Q_k) \phi'_j(\xi^Q_k) \frac{1}{J(\xi^Q_k)} \\
   &+ \phi_i^\prime (\xi^Q_k) \underline{\underline{\mathcal{G}}}^\mathrm{D1}(\xi^Q_k) \phi_j(\xi^Q_k) \\
   & +\phi_i(\xi^Q_k) \underline{\underline{D}}^\mathrm{D2}(\xi^Q_k) \phi^\prime_j(\xi^Q_k) \\
   & \phi_i(\xi^Q_k) \left[\underline{\underline{\mathcal{G}}}^\mathrm{I}(\xi^Q_k) 
   + \underline{\underline{\mathcal{G}}}^\mathrm{D2}(\xi^Q_k) \right] \phi_j(\xi^Q_k) J(\xi^Q_k)
   \Big\}w^Q_k \\
   %
   \underline{\underline{K}}_{ij} =
   \sum_{k=1}^{n^Q} 
   \Big\{& 
   \phi'_i(\xi^Q_k) \left[
   \underline{\underline{C}}(\xi^Q_k) 
   + \underline{\underline{D}}^\mathrm{D1}(\xi^Q_k) 
    \right] \phi'_j(\xi^Q_k) \frac{1}{J(\xi^Q_k)} \\
   & + \phi'_i(\xi^Q_k) \left[
   \underline{\underline{K}}^\mathrm{E1}(\xi^Q_k) 
   + \underline{\underline{K}}^\mathrm{D1}(\xi^Q_k) 
    \right] \phi_j(\xi^Q_k) \\
   & \phi_i(\xi^Q_k) \left[ 
   \underline{\underline{\mathcal{P}}}^\mathrm{E2}(\xi^Q_k) +
   \underline{\underline{\mathcal{P}}}^\mathrm{D2}(\xi^Q_k) \right]
   \phi'_j(\xi^Q_k) \\
   & + \phi_i(\xi^Q_k) \left[
   \underline{\underline{\mathcal{K}}}^\mathrm{E2}(\xi^Q_k)
   +\underline{\underline{\mathcal{K}}}^\mathrm{D2}(\xi^Q_k) 
   \right]\phi_j(\xi^Q_k) J(\xi^Q_k)
   \Big\} w^Q_k 
   \end{aligned}

for all :math:`i,j \in\{1,2, \ldots, P\}` and
:math:`\underline{\underline{M}}_{ij},\underline{\underline{G}}_{ij},\underline{\underline{K}}_{ij} \in \mathbb{R}^{6 \times 6}`.
These matrices define the full matrices for a single beam LSFE. For
example, the mass matrix is assembled as

.. math::

   \underline{\underline{M}} =
   \begin{bmatrix}
   \underline{\underline{M}}_{11}&
   \underline{\underline{M}}_{12}& \ldots &
   \underline{\underline{M}}_{1P}\\
   \underline{\underline{M}}_{21} &
   \underline{\underline{M}}_{22}&
   \ldots &
   \underline{\underline{M}}_{2P}\\
   \vdots & \vdots & \vdots & \vdots \\
   \underline{\underline{M}}_{P1}&
   \underline{\underline{M}}_{P2}& \ldots &
   \underline{\underline{M}}_{PP}\\
   \end{bmatrix}

where
:math:`\underline{\underline{M}} \in \mathbb{R}^{6 P \times 6 P}`;
similarly for
:math:`\underline{\underline{G}}, \underline{\underline{K}} \in \mathbb{R}^{6 P \times 6 P}`.

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Deville-etal:2002

      Deville, M. O., P. F. Fischer, and E. H. Mund. 2002. *High-Order
      Methods for Incompressible Fluid Flow*. Cambridge University
      Press.
