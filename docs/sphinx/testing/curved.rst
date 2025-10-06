Static bending of a curved beam
-------------------------------

In this benchmark problem, we examine the tip deflection of a tip-loaded curved beam, which is described in detail in [@Bathe-Bolourchi:1979].
The beam lies in the x-z plane with a tip force of 600 lbs acting in the 
positive y direction.  Beam is curved on a radius of 100 inches in the positive 
x direction through a 45 degree arc.

Geometric properties

- Cross section: 1 in x 1 in

Material properties

- Elastic Modulus, E = 10**7 psi
- Poisson Ratio, nu = 0

Diagonal stiffness matrix entries

- kGA = 41.6667E+5
- kGA = 41.6667E+5
- EA = 100.E+5
- EI_1 = 8.333E+5
- EI_2 = 8.333E+5
- GJ = 8.333E+5

Results

Tip displacement reported by [@Bathe-Bolourchi:1979]:
(-13.4, 53.4, -23.5)


The following table show the tip displacement calculated with a single Legendre Spectral Finite Element (LSFE) with the number of nodes shown. N-Point Gauss-Legendre quadrature was employed, where N is the number of nodes. The results show the expected spectral convergence to a converged solution.

+----------------------+-------------------------+-------------------------+-------------------------+
| Number of Nodes (N)  |           u1            |           u2            |           u3            |
+======================+=========================+=========================+=========================+
|          3           | -1.390261363798350E+00  |   2.264633130098760E+01 |   3.168614893736530E+01 |
+----------------------+-------------------------+-------------------------+-------------------------+
|          6           | -2.362791403924360E+01  |   1.372486827629740E+01 |   5.356617784705440E+01 |
+----------------------+-------------------------+-------------------------+-------------------------+
|          9           | -2.364597006589220E+01  |   1.372080411971940E+01 |   5.357763167143680E+01 |
+----------------------+-------------------------+-------------------------+-------------------------+
|         12           | -2.364597159890480E+01  |   1.372080060483760E+01 |   5.357764506053460E+01 |
+----------------------+-------------------------+-------------------------+-------------------------+
|         15           | -2.364597160090470E+01  |   1.372080060400580E+01 |   5.357764506269870E+01 |
+----------------------+-------------------------+-------------------------+-------------------------+
|         18           | -2.364597160090480E+01  |   1.372080060400580E+01 |   5.357764506269890E+01 |
+----------------------+-------------------------+-------------------------+-------------------------+


_References_

   .. container:: csl-entry
      :name: ref-Bathe-Bolourchi:1979

      Bathe, K.-J., and S. Bolourchi. 1979. “Large displacement analysis
      of three-dimensional beam structures."
      *International Journal for Numerical Methods in Engineering* 14:
      961-986.

