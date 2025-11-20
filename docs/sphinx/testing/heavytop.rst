.. _sec-heavytop-verification:

Heavy Top Dynamic Problem
-------------------------

For this benchmark problem, we use the heavy top to evaluate the convergence rate of the Kynema
Lie-group generalized-alpha time integrator. The problem consists of a rotating body, fixed to the
ground by a spherical joint, as detailed in [@Bruls-etal:2012]. The formulation of the heavy top
problem as implemented in Kynema is described in :ref:`sec-heavy-top`.

Our benchmark solution is computed using a Kynema simulation with a time step size of :math:`10^{-6}`
seconds; we examine the vertical displacement :math:`u_3` at the end of two seconds. The following
figure shows the convergence of the relative error in the vertical displacement as a function of
time step size, where we observe the expected second-order convergence of the method.

.. figure:: images/heavy_top_u3_convergence_plot.png
   :align: center
   :width: 80%

   Second-order convergence of the relative error in vertical displacement :math:`u_3` at :math:`t = 2.0` s
   for the heavy top problem. The benchmark solution is computed with :math:`\Delta t = 10^{-6}` s.

.. note::

   This benchmark is included as a regression test in the Kynema code base.
   The test is implemented in:

     https://github.com/kynema/kynema/blob/main/tests/regression_tests/regression/verification/dynamic_heavy_top.cpp

.. container:: csl-entry
   :name: ref-Bruls-etal:2012

   Brüls, O., A. Cardona, and M. Arnold. 2012. “Lie Group Generalized-:math:`\alpha` Time Integration
   For Constrained Flexible Multibody Systems.” *Mechanism and Machine Theory* 48: 121–137.
