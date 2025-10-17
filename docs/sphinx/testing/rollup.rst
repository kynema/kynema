.. _sec-rollup:

Static rollup of an isotropic beam
----------------------------------

In this benchmark test, we examine the static pure bending of an isotropic beam of 10 units in
length. Under the application of a specific moment, the beam rolls up into a circular curve and
the exact analytical tip displacement (in the direction of the beam length) is -10. Details of
the benchmark problem are described in [@Simo-VuQuoc:1986]. The following table shows the tip
displacement calculated with a single Legendre Spectral Finite Element (LSFE) with the number
of nodes shown. N-Point Gauss-Legendre quadrature was employed, where N is the number of nodes.
The results show the expected spectral convergence to the exact solution.

+------------+--------------------+
| # of nodes | Tip Displacement   |
+============+====================+
| 3          | -15.6418448333196  |
+------------+--------------------+
| 6          | -10.9854783002967  |
+------------+--------------------+
| 9          | -10.058403575783   |
+------------+--------------------+
| 12         | -10.0000013123339  |
+------------+--------------------+
| 15         | -10.0000000000645  |
+------------+--------------------+


.. note::

   This benchmark is included as a regression test in the Kynema code base.
   The test is implemented in:

     https://github.com/kynema/kynema/blob/main/tests/regression_tests/regression/verification/static_isotropic_beam_rollup.cpp


.. container:: csl-entry
   :name: ref-Simo-VuQuoc:1985


   Simo, J. C., and L. Vu Quoc. 1986. “A three-dimensional
   finite-strain rod model. Part II."
   *Computer Methods in Applied Mechanics and Engineering* 58:
   79-116.


