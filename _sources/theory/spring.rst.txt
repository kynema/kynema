.. _sec-spring:

Geometrically nonlinear spring
------------------------------

In this section we describe the variables and unconstrained governing
equations for a two-node spring that is constitutively linear and
geometrically nonlinear. The spring input properties are the spring constant :math:`k` and its unstretched length :math:`L`. The spring
element has reference position defined by its endpoint positions,

.. math::

   \begin{aligned}
    \underline{q}^\mathrm{r} = 
   \begin{bmatrix}
     \underline{x}_1^\mathrm{r} \\
     \underline{x}_2^\mathrm{r} 
   \end{bmatrix}
   \end{aligned}

where :math:`\underline{q}^\mathrm{r} \in \mathbb{R}^6`,
and displacement is denoted

.. math::

   \begin{aligned}
    \underline{q} = 
   \begin{bmatrix}
     \underline{u}_1 \\
     \underline{u}_2 
   \end{bmatrix}
   \end{aligned}

The spring internal force contribution to Eq. :eq:`residual1` is

.. math::

   \underline{g} = \begin{bmatrix}
   \underline{f}^\mathrm{sp} \\
   -\underline{f}^\mathrm{sp}
   \end{bmatrix}

where :math:`\underline{g} \in \mathbb{R}^6` and

.. math:: \underline{f}^\mathrm{sp} = -k \frac{\underline{r} }{| \underline{r} |} \left( | \underline{r} | - L \right)

with 

.. math:: \underline{r} = \underline{x}_2^\mathrm{r} + \underline{u}_2 - \underline{x}_1^\mathrm{r}  - \underline{u}_1

Variation of the force equation provides the stiffness contribution
to the generalized-:math:`\alpha` iteration matrix (See Eq. :eq:`iteration`):

.. math::

   \underline{\underline{K}} =  \begin{bmatrix}
   \underline{\underline{A}} & -\underline{\underline{A}} \\
   - \underline{\underline{A}} & \underline{\underline{A}}
   \end{bmatrix}

and

.. math::

   \underline{\underline{A}} =  k \left( \frac{L}{|\underline{r} |} - 1\right) \underline{\underline{I}}
   - \frac{k L}{|\underline{r}|^3}\left( \widetilde{r} \widetilde{r} + \underline{\underline{I}} (\underline{r}^T \underline{r} ) \right)
