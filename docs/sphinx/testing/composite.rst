.. _sec-composite:

Static bending of a straight composite beam
-------------------------------------------

In this benchmark problem we examine the tip deflection of a composite beam
under a tip load.  The problem is described in detail in [@Wang-etal:2017].
The problem was solved with a single Legendre Spectral Finite Element (LSFE)
with N-point Gauss-Legendre quadrature where N is the number of nodes.  The
tip-displacement results, shown in the table below, are the same (to within
machine precision) of those produced with BeamDyn.


+---------+----------------------+---------------------+--------------------+
| # nodes |       u1             |     u2              |    u3              |
+=========+======================+=====================+====================+
| 3       | -5.04479565690760E-02|-3.91038892649929E-02|9.11237052941376E-01|
+---------+----------------------+---------------------+--------------------+
| 6       | -9.02727162973492E-02|-6.47489265976265E-02|1.22973611590669E+00|
+---------+----------------------+---------------------+--------------------+
| 9       | -9.02726627568822E-02|-6.47488486249749E-02|1.22973648291939E+00|
+---------+----------------------+---------------------+--------------------+
| 12      | -9.02726627566302E-02|-6.47488486259037E-02|1.22973648292371E+00|
+---------+----------------------+---------------------+--------------------+
| 15      | -9.02726627566299E-02|-6.47488486259036E-02|1.22973648292371E+00|
+---------+----------------------+---------------------+--------------------+
| 18      | -9.02726627566296E-02|-6.47488486259039E-02|1.22973648292371E+00|
+---------+----------------------+---------------------+--------------------+


.. note::

   This benchmark is included as a regression test in the Kynema code base.
   The test is implemented in:

     https://github.com/kynema/kynema/blob/main/tests/regression_tests/regression/verification/static_composite_beam_bending.cpp


.. container:: csl-entry
   :name: ref-Wang-etal:2017

   Wang, Q., M.A. Sprague, J. Jonkman, N. Johnson, B. Johnkman. 2017. “BeamDyn:
   a high-fidelity wind turbine blade solver in the FAST modular framework."
   *Wind Energy* 20: 1439-1462.
   https://onlinelibrary.wiley.com/doi/pdf/10.1002/we.2101

