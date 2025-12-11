# NREL5MW Turbine

This test shows how to run a simulation of a full wind turbine model based on a realistic specification of the IEA15MW turbine.
It uses Kynema's high level "Turbine" interface for easy definition of this common problem configuration.

To build this test, first ensure that Kynema with ROSCO has been installed somewhere that is discoverable by CMake (if using the Spack package manager, run `spack load kynema`).
Next, create a build directory and from there run `cmake ../` and `make`.
When running this test, you will need to be in the same directory as the included YAML input file, the DISCON file, and the associated Cp_Ct_Cq file.
