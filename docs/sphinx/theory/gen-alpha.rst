.. _sec-gen-alpha:

Time integration: Generalized-Alpha
-----------------------------------

Kynema time-integration is achieved with a Lie-group
generalized-:math:`\alpha` algorithm, described in detail by
[@Bruls-etal:2012]. This algorithm is second-order time accurate and
summarized here for completeness. For a system represented as
Eq. :eq:`residual` and with linearized form shown in
Eq. :eq:`variation`, the generalized-:math:`\alpha`
algorithm can be written as in
Algorithm 1, where an
:math:`n` superscript denotes evaluation at the :math:`n`\ th time step,
and :math:`\rho_\infty \in[0,1]` is the numerical damping with
:math:`\rho_\infty = 1` being no damping and :math:`\rho_\infty=0` being
maximum damping. The term :math:`\underline{a}^{n} \in \mathbb{R}^k` is
an auxiliary variable [@Arnold-Bruls:2007] for computation and is
initialized as :math:`\underline{a}^{0}= \dot{\underline{v}}^0`. The
algorithm includes the left and right preconditioning matrices described
in [@Bottasso-etal:2008],

.. math::

   \begin{aligned}
   \underline{\underline{D}}_L\left( \beta \Delta t^2\right) =
   \begin{bmatrix}
   \beta \Delta t^2 \underline{\underline{I}}_k & \underline{\underline{0}} \\
   \underline{\underline{0}} & \underline{\underline{I}}_m
   \end{bmatrix}
   \quad
   \underline{\underline{D}}_R\left( \beta \Delta t^2\right) =
   \begin{bmatrix}
   \underline{\underline{I}}_k & \underline{\underline{0}} \\
   \underline{\underline{0}} & \frac{1}{\beta \Delta t^2} \underline{\underline{I}}_m
   \end{bmatrix}
   \end{aligned}

respectively, where
:math:`\underline{\underline{D}}_L,\underline{\underline{D}}_R \in \mathbb{R}^{k\times k}`,
:math:`\underline{\underline{I}}_p` is the :math:`p \times p` identity
matrix, and :math:`\beta` is a parameter defined in Algorithm 1.

Important to the time-integration scheme is the exponential map, which
can be written for :math:`k_n` nodes

.. math::

   \exp_{\mathbb{R}^3\times\mathrm{SO(3)}} 
   \begin{bmatrix} 
   \underline{u}_1 \\
   \underline{\psi}_1 \\
   \vdots \\
   \underline{u}_{k_n} \\
   \underline{\psi}_{k_n} 
   \end{bmatrix} 
   = \begin{bmatrix} 
   \underline{u}_1 \\ 
   \exp_\mathrm{SO(3)} \left(\widetilde{\psi}_1\right) \\
   \vdots \\
   \underline{u}_{k_n} \\ 
   \exp_\mathrm{SO(3)} \left(\widetilde{\psi}_{k_n}\right) \\
   \end{bmatrix}

where

.. math:: \exp_\mathrm{SO(3)}\left(\widetilde{\psi}\right)  = \underline{\underline{I}} + \frac{\sin \phi}{\phi} \widetilde{\psi} + \frac{1-\cos \phi}{\phi^2} \widetilde{\psi}\widetilde{\psi}

:math:`\underline{\underline{I}} \in \mathbb{R}^{3 \times 3}` is the
identity matrix, :math:`\underline{\psi}` is the Cartesian rotation
vector, :math:`\phi = || \underline{\psi} ||`, and the tilde operator is
defined for a vector :math:`\underline{a}\in \mathbb{R}^3` with entries
:math:`a_i` as

.. math::

   \widetilde{a} =
   \begin{bmatrix}
   0 & -a_3 & a_2 \\
   a_3 & 0 & -a_1 \\
   -a_2 & a_1 & 0
   \end{bmatrix}



----------------------------------------------

**Algorithm 1:** Lie group generalized-α time integration for one time step [@Bruls-etal:2012]

----------------------------------------------

**Require:** :math:`\underline{q}^n`, :math:`\underline{v}^n`, :math:`\dot{\underline{v}}^n`, :math:`\underline{a}^n`, :math:`t^n`, :math:`\Delta t`, :math:`\rho_\infty`, atol, rtol, :math:`i_\mathrm{max}`, :math:`k`, :math:`m`

      :math:`\alpha_m = \frac{2 \rho_\infty - 1}{\rho_\infty+1}`

      :math:`\alpha_f = \frac{\rho_\infty}{\rho_\infty+1}`

      :math:`\gamma = 0.5 + \alpha_f - \alpha_m`

      :math:`\beta = 0.25 \left( \gamma + 0.5\right)^2`

      :math:`\beta^\prime = \frac{1-\alpha_m}{\Delta t^2 \beta (1-\alpha_f)}`

      :math:`\gamma^\prime = \frac{\gamma}{\Delta t \beta}`

      :math:`t^{n+1} := t^n + \Delta t`

      :math:`\dot{\underline{v}}^{n+1} := \underline{0}`

      :math:`\underline{\lambda}^{n+1} := \underline{0}`

      :math:`\underline{a}^{n+1} := (\alpha_f \dot{\underline{v}}^{n} - \alpha_m \underline{a}^n)/(1-\alpha_m)`

      :math:`\underline{v}^{n+1} := \underline{v}^n + \Delta t (1-\gamma) \underline{a}^n + \gamma \Delta t \underline{a}^{n+1}`

      :math:`\Delta \underline{q}^n := \underline{v}^n+(0.5-\beta) \Delta t \underline{a}^n + \beta \Delta t \underline{a}^{n+1}`

      :math:`\mathrm{err} := 2.0`    *Initialize with any* :math:`\mathrm{err} > 1.0` 

      :math:`i := 0` 

      While :math:`i \le i_\mathrm{max}` and :math:`\mathrm{err} > 1.0` do

          :math:`i := i+1`

          :math:`\underline{q}^{n+1} := \exp_{\underline{\underline{R}}^3\times \mathrm{SO(3)}} ( \Delta t \Delta \underline{q}^n ) \circ \underline{q}^n`

          Calculate :math:`\underline{r}^{n+1}`; See Eq. :eq:`residual1` 

          Calculate :math:`\underline{\underline{S}}_t^{n+1}`; See Eq. :eq:`iteration`  

          Solve :math:`\underline{\underline{D}}_L(\beta \Delta t^2) \underline{\underline{S}}_t^{n+1} \underline{\underline{D}}_R(\beta \Delta t^2) \begin{bmatrix} \Delta \underline{x} \Delta \underline{\lambda} \end{bmatrix}= -\underline{\underline{D}}_L(\beta \Delta t^2) \underline{r}^{n+1}`

          :math:`\begin{bmatrix} \Delta \underline{x} \Delta \underline{\lambda} \end{bmatrix} := \underline{\underline{D}}_R(\beta \Delta t^2) \begin{bmatrix} \Delta \underline{x} \Delta \underline{\lambda} \end{bmatrix}`

          :math:`\Delta \underline{q}^n := \Delta \underline{q}^n + \Delta \underline{x}/\Delta t^n`

          :math:`\underline{v}^{n+1} := \underline{v}^{n+1} + \gamma^\prime \Delta \underline{x}`

          :math:`\dot{\underline{v}}^{n+1} := \dot{\underline{v}}^{n+1} + \beta^\prime \Delta \underline{x}`

          :math:`\underline{\lambda}^{n+1} := \underline{\lambda}^{n+1} + \Delta \underline{\lambda}`

          :math:`\mathrm{err} := \sqrt{ \frac{1}{k + m} \left( \sum_{i=1}^{k} \left( \frac{ \Delta x_i }{  \mathrm{atol} + \mathrm{rtol} \left| \Delta t \Delta q_i^n \right| } \right)^2 + \sum_{i=1}^{m} \left( \frac{ \Delta \lambda_i }{  \mathrm{atol} + \mathrm{rtol} \left| \lambda_i^{n + 1} \right| } \right)^2 \right) }` See [@Arnold-Hante:2017]

      **end while**
    
      :math:`\underline{a}^{n+1} := \underline{a}^{n+1} + \dot{\underline{v}}^{n+1}\left( 1 - \alpha_f\right) / \left( 1 - \alpha_m\right)`

      **Return:** :math:`t^{n+1}`, :math:`\underline{q}^{n+1}`,
      :math:`\underline{v}^{n+1}`, :math:`\dot{\underline{v}}^{n+1}`,
      :math:`\underline{\lambda}^{n+1}`, :math:`\underline{a}^{n+1}`

--------------------------------------------------------------

The so-called iteration matrix,
:math:`\underline{\underline{S}}_t \in \mathbb{R}^{(k+m)\times (k+m)}`,
for Eqs. :eq:`residual` and :eq:`variation` can be written

.. math:: \underline{\underline{S}}_t = \begin{bmatrix}
   \underline{\underline{M}} \beta'+\underline{\underline{G}} \gamma' + \left(\underline{\underline{K}} + \underline{\underline{K}}^\Phi\right)\, \underline{\underline{T}}_{\mathbb{R}^3\times \mathrm{SO(3)}}^T(\Delta t \Delta q) & \underline{\underline{B}}^T \\
   \underline{\underline{B}}\,\underline{\underline{T}}_{\mathbb{R}^3\times \mathrm{SO(3)}}^T(\Delta t \Delta q)                     & \underline{\underline{0}}
   \end{bmatrix}
   :label: iteration

where the tangent matrix,
:math:`\underline{\underline{T}}_{\mathbb{R}^3\times \mathrm{SO(3)}}(\underline{\psi}) \in \mathbb{R}^{k\times k}`,
is

.. math::

   \begin{aligned}
   \underline{\underline{T}}_{\mathbb{R}^3\times \mathrm{SO(3)}}(\underline{\psi})  = 
   \begin{bmatrix} 
   \underline{\underline{I}} & \underline{\underline{0}}                            &        & \cdots &  \underline{\underline{0}}\\ 
   \underline{\underline{0}} & \underline{\underline{T}}_{\mathrm{SO(3)}}(\underline{\psi}_1) &        &        & \\
          &                                   &\ddots  &        & \\
          &                                   &        & \underline{\underline{I}} & \underline{\underline{0}} \\
    \underline{\underline{0}}&   \cdots                          &        & \underline{\underline{0}} & \underline{\underline{T}}_{\mathrm{SO(3)}}(\underline{\psi}_{k_n}) 
   \end{bmatrix}
   \end{aligned}

and the variation of the virtual rotation is related to the variation of
the Cartesian rotation vector as

.. math::

   \begin{aligned}
   \delta \underline{\theta} = \underline{\underline{T}}_{\mathrm{SO(3)}}^T(\underline{\psi}) \delta \underline{\psi}
   \end{aligned}

with [@Geradin-Cardona:2001]

.. math::

   \begin{aligned}
   \underline{\underline{T}}_{\mathrm{SO(3)}}(\underline{\psi}) = \underline{\underline{I}} 
   + \left(\frac{\cos ||\underline{\psi}|| -1}{||\underline{\psi}||^2} \right) \widetilde{\psi}
   +\left(1- \frac{\sin ||\underline{\psi}||}{||\underline{\psi}||}\right) 
   \frac{\widetilde{\psi} \widetilde{\psi}}{||\underline{\psi}||^2}
   \end{aligned}

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Arnold-Bruls:2007

      Arnold, M., and O. Brüls. 2007. “Convergence of the
      Generalized-:math:`\alpha` Scheme for Constrained Mechanical
      Systems.” *Multibody System Dynamics* 18: 185–202.

   .. container:: csl-entry
      :name: ref-Arnold-Hante:2017

      Arnold, M., and S. Hante. 2017. “Implementation Details of a
      Generalized-:math:`\alpha` Differential-Algebraic Equation Lie
      Group Method.” *Journal of Computational and Nonlinear Dynamics*
      2: 021002.

   .. container:: csl-entry
      :name: ref-Bottasso-etal:2008

      Bottasso, C. L., D. Dopico, and L. Trainelli. 2008. “On the
      Optimal Scaling of Index Three DAEs in Multibody Dynamics.”
      *Multibody System Dynamics* 19: 3–20.

   .. container:: csl-entry
      :name: ref-Bruls-etal:2012

      Brüls, O., A. Cardona, and M. Arnold. 2012. “Lie Group
      Generalized-:math:`\alpha` Time Integration Fo Constrained
      Flexible Multibody Systems.” *Mechanism and Machine Theory* 48:
      121–137.

   .. container:: csl-entry
      :name: ref-Geradin-Cardona:2001

      Géradin, M., and A. Cardona. 2001. *Flexible Multibody Dynamics: A
      Finite Element Approach*. Chichester: John Wiley & Sons.
