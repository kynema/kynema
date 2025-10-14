.. _sec-twisted:

Static bending of a twisted beam
--------------------------------

For this benchmark problem, we examine the static tip deflection of a straight, but twisted beam. Details are provided in
[@Wang-etal:2017].  The benchmark solution is a highly refined solid-element model solved in ANSYS.  The following table show the tip displacement calculated with a single Legendre Spectral Finite Element (LSFE) with the number of nodes shown. N-Point Gauss-Legendre quadrature was employed, where N is the number of nodes. The results show the expected spectral convergence a solution that is in good agreement with the ANSYS solid-element solution.  The results are also identical, within machine precision, to those produced with BeamDyn.


ANSYS Solid-Element Solution: u1=-1.134192, u2= -1.714467, u3=-3.58423


+----------------------+------------------------+------------------------+------------------------+
| Number of Nodes (N)  |          u1            |          u2            |          u3            |
+======================+========================+========================+========================+
|          3           | -1.142067384519760E-01 | -1.984490462570690E-01 | -1.296500991573600E+00 |
+----------------------+------------------------+------------------------+------------------------+
|          6           | -1.124292284062730E+00 | -1.690365192280430E+00 | -3.573215914495870E+00 |
+----------------------+------------------------+------------------------+------------------------+
|          9           | -1.141513707963730E+00 | -1.718037758175690E+00 | -3.593187973051160E+00 |
+----------------------+------------------------+------------------------+------------------------+
|         12           | -1.141525993948990E+00 | -1.718053339951670E+00 | -3.593213227112630E+00 |
+----------------------+------------------------+------------------------+------------------------+
|         15           | -1.141525999539420E+00 | -1.718053346463880E+00 | -3.593213240269520E+00 |
+----------------------+------------------------+------------------------+------------------------+
|         18           | -1.141525999540470E+00 | -1.718053346464930E+00 | -3.593213240271990E+00 |
+----------------------+------------------------+------------------------+------------------------+


   .. container:: csl-entry
      :name: ref-Wang-etal:2017

      Wang, Q., M. A. Sprague, J. Jonkman, N. Johnson, and B. Jonkman.
      2017. “BeamDyn: A High-Fidelity Wind Turbine Blade Solver in the
      FAST Modular Framework.” *Wind Energy* 20: 1439-1462.
      https://onlinelibrary.wiley.com/doi/pdf/10.1002/we.2101
