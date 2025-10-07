.. _sec-quadrature:

Notes on quadrature
^^^^^^^^^^^^^^^^^^^

In Kynema, we represent beams that can have highly variable properties
along the length and we use a single high-order element for the whole
beam. Material properties, including the sectional mass matrix
:math:`\underline{\underline{M}}^*`, stiffness matrix
:math:`\underline{\underline{C}}^*`, and twist :math:`\tau`, are defined
by the user at stations located along the beam reference line. In typical
finite-element beam implementations, :math:`P`-point Gauss-Legendre
quadrature is commonly used for evaluating the underlying weak-form
integrals. While that is often sufficient for uniform or linearly varying
properties, it can be inadequate for highly variable material properties.
Kynema provides users the option of Gauss-Legendre (GL) or
Gauss-Lobatto-Legendre (GLL) quadrature.  GL or GLL quadrature can be
applied over the whole beam with a user-specified number of points, or in
a composite fashion between user-defined material points.  In the latter,
a composite GLL quadrature with 2 points would be equivalent to using
trapezoid rule integration over the full blade with an integration point
at each user-defined material station.  
