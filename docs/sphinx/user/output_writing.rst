NetCDF Output Writing
=====================

This guide explains how Kynema writes simulation results to NetCDF files, the data schemas used
for node state and time-series outputs, and how to inspect and work with the generated files.
NetCDF (Network Common Data Form) is a self-describing, machine-independent binary format
widely used in scientific computing for storing array-oriented data with rich metadata.

Overview
--------

Kynema automatically generates NetCDF output files when using the high-level
interface APIs provided in ``src/interfaces``. The output system provides two
complementary (however independent) data streams:

1. **Node State Data**: Full structural state snapshots for all nodes at each timestep,
   including position, displacement, velocity, acceleration, and forces. Written via
   ``NodeStateWriter`` in ``src/utilities/netcdf/node_state_writer.hpp``.

2. **Time-Series Data**: Custom scalar or vector quantities tracked over time, such as
   controller outputs, constraint forces, or user-defined metrics. Written via
   ``TimeSeriesWriter`` in ``src/utilities/netcdf/time_series_writer.hpp``.

These files are the foundation for post-processing workflows, including the VTK visualization
pipeline described in :doc:`vtk_visualization`.

Prerequisites
-------------

Optionally install the NetCDF command-line utilities so you can inspect generated files quickly. The
``ncdump`` tool (packaged as ``netcdf-bin`` on Debian/Ubuntu) prints headers and variable data in
plain text:

.. code-block:: bash

    sudo apt-get install netcdf-bin

No additional runtime steps are required; Kynema links the NetCDF-C++ library through
``cmake/Dependencies.cmake`` when you enable output writing in a build.

Output Files and Naming Conventions
------------------------------------

Default File Naming
~~~~~~~~~~~~~~~~~~~

When simulations are run via Kynema's interface APIs, output files are automatically generated
in the working directory or a specified output path. The naming convention depends on the
interface used:

- **TurbineInterface**: ``turbine_interface.nc`` (node state) + optional time-series file
- **CFDInterface**: ``cfd_interface.nc`` (node state) + optional time-series file  
- **BladeInterface**: ``blade_interface.nc`` (node state) + optional time-series file

Customizing Output Paths
~~~~~~~~~~~~~~~~~~~~~~~~~

You can customize the output file path using the builder pattern. For example, in the
``BladeInterface``:

.. code-block:: cpp

    auto interface = TurbineInterfaceBuilder{}
        .Outputs().SetOutputFilePath("results/run01")
        .Build();

CFDInterface exposes ``SetOutputFile(...)`` with the same semantics.

Node State Data Schema
----------------------

The ``NodeStateWriter`` class manages structured storage of nodal quantities at each timestep.
The schema uses an unlimited time dimension, allowing continuous writing as the simulation
progresses.

Dimensions
~~~~~~~~~~

- **time**: Unlimited dimension representing simulation timesteps
- **nodes**: Fixed dimension equal to the number of nodes in the model

State Components
~~~~~~~~~~~~~~~~

Each prefix groups a set of variables with dimensions ``[time, nodes]``:

- ``x_*`` — global position (x, y, z) and orientation quaternion components (w, i, j, k)
- ``u_*`` — displacement from the reference configuration - translation (x, y, z)  and rotation (w, i, j, k)
- ``v_*`` — translational (x, y, z) and rotational (i, j, k) velocity
- ``a_*`` — translational (x, y, z) and rotational (i, j, k) acceleration
- ``f_*`` — nodal forces (x, y, z) and moments (i, j, k)

Writing Options
~~~~~~~~~~~~~~~

``NodeStateWriter`` accepts optional tuning parameters:

- ``enabled_state_prefixes`` filters the prefixes written to disk (default is to write all components
  ``{"x", "u", "v", "a", "f"}``).
- ``buffer_size`` batches timesteps before a flush; set to a positive integer for higher throughput or to ``0`` for immediate writes.
- Chunk sizes are derived from the buffer setting, so most workloads do not need manual tuning.

.. code-block:: cpp

    // Write only displacement and velocity at 50-timestep intervals
    util::NodeStateWriter writer(
        "output.nc",
        true,                      // create new file
        num_nodes,
        {"u", "v"},                // enabled components
        50                         // buffer 50 timesteps before auto-flush
    );

Time-Series Data Schema
------------------------

``TimeSeriesWriter`` captures scalar or vector channels indexed by an unlimited ``time``
dimension. The first write defines the corresponding ``{name}_dimension`` and the variable
shape, so the writer naturally supports custom schemas without a separate setup step. When a
time-series file is configured, ``outputs.GetTimeSeriesWriter()`` returns a live writer;
otherwise it is null.

Turbine-oriented interfaces pre-build a catalog via ``TurbineTimeSeriesSchema`` (see
``time-series/src/interfaces/turbine/turbine_time_series_schema.hpp``) to guarantee consistent
channel ordering for aero, structural, and controller quantities. Additional channels can be
registered ad hoc:

.. code-block:: cpp

    // Called from the simulation loop
    outputs.WriteValueAtTimestep(timestep, "generator_torque", torque);
    if (auto& series = outputs.GetTimeSeriesWriter()) {
        series->WriteValuesAtTimestep(
            "blade_pitch_angles", timestep, std::array{pitch_a, pitch_b, pitch_c}
        );
    }

Turbine Time-Series Channels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Turbine-oriented interface ``TurbineInterface`` publishes a standardized set of channels for post-processing. The
schema is dynamically constructed to match the configured turbine and available components:

- Basic simulation:
    - Time (s): simulation time in seconds.
    - ConvIter (-): number of convergence iterations performed for the timestep.
    - ConvError (-): final convergence error metric for the timestep.
    - Azimuth (deg): rotor azimuth angle.
    - RotSpeed (rpm): rotor rotational speed.
    - YawPzn (deg): nacelle yaw position angle.

- Tower top/base:
    - YawBrTDxt/TDyt/TDzt (m): tower-top displacement components [x,y,z] in the global frame.
    - YawBrTVxp/TVyp/TVzp (m_s): tower-top translational velocity components [x,y,z] in the global frame.
    - YawBrTAxp/TAyp/TAzp (m_s^2): tower-top translational acceleration components [x,y,z] in the global frame.

    - TwrBsFxt/Fyt/Fzt (kN): tower-base reaction forces [Fx,Fy,Fz] in the global frame.
    - TwrBsMxt/Myt/Mzt (kN-m): tower-base reaction moments [Mx,My,Mz] in the global frame.

- Thrust
    - RotThrust (kN): axial rotor thrust.

- Per-blade group:
    - RootFxr/Fyr/Fzr (N): blade root forces [Fx,Fy,Fz] in the rotor frame.
    - RootMxr/Myr/Mzr (N-m): blade root moments [Mx,My,Mz] in the rotor frame.
    - BldPitch (deg): blade pitch angle per blade.
    - TipTVXg/TVYg/TVZg (m_s): blade tip translational velocity components [X,Y,Z] in the global frame.
    - TipRVXg/RVYg/RVZg (deg_s): blade tip rotational velocity components [X,Y,Z] in the global frame.
    - Per-blade stride: 13 channels per blade within this group.

- Controller (if enabled):
    - GenTq (kN-m): generator torque (controller output, if enabled).
    - GenPwr (kW): generator electrical power (controller output, if enabled).

- Hub inflow:
    - WindHubVelX/Y/Z (m_s): hub-height inflow wind velocity components [X,Y,Z] in the global frame.

- Aerodynamics (per body and section):
    - Vrel (m_s): section relative flow speed.
    - Alpha (deg): section angle of attack.
    - Cn/Ct/Cm (-): section aerodynamic coefficients (normal, tangential, pitching moment).
    - Fxi/Fyi/Fzi (N_m): section internal force terms in the local element frame.
    - Mxi/Myi/Mzi (N_m): section internal moment terms in the local element frame.
    - Labels use AB{blade}N{node} with a 3‑digit node label (e.g., AB1N001).
    - Per-section stride: 11 channels per section within this group.

Note: Channel presence and counts are configuration-dependent (e.g., number of blades, aero bodies,
controller availability). You can append custom channels using the time-series writer during the
simulation loop.

Inspecting NetCDF Files
------------------------

Using ncdump
~~~~~~~~~~~~

The ``ncdump`` utility provides a quick way to inspect NetCDF file structure and contents
from the command line.

**View file metadata and structure:**

.. code-block:: bash

    ncdump -h output.nc

This displays dimensions, variables, and attributes without dumping the full data arrays.

These variables are the inputs consumed by the VTK conversion workflow.

**Example output for node state file:**

.. code-block:: text

    netcdf output {
    dimensions:
        time = UNLIMITED ; // (1001 currently)
        nodes = 51 ;
    variables:
        double x_x(time, nodes) ;
        double x_y(time, nodes) ;
        double x_z(time, nodes) ;
        double x_w(time, nodes) ;
        double x_i(time, nodes) ;
        double x_j(time, nodes) ;
        double x_k(time, nodes) ;
        double u_x(time, nodes) ;
        // ... additional variables
    }

**Extract specific variable data:**

.. code-block:: bash

    # Dump all data for displacement in x-direction
    ncdump -v u_x output.nc

    # Dump data for multiple variables
    ncdump -v u_x,u_y,u_z output.nc

Examples from Regression Tests
------------------------------

Several regression tests exercise NetCDF writing and provide reference configurations:

- ``NetCDFOutputsWriterTest.SpringMassSystemOutputs`` in
  ``tests/regression_tests/regression/test_netcdf_outputs_writer.cpp`` shows direct use of
  ``NodeStateWriter`` for a minimal 2-node spring–mass model.
- ``CFDInterfaceTest.FloatingPlatform`` in
  ``tests/regression_tests/interfaces/test_cfd_interface.cpp`` enables interface-driven output
  and validates displacement values read back from ``cfd_interface.nc``.
- ``TurbineInterfaceTest.IEA15_Structure`` in
  ``tests/regression_tests/interfaces/test_turbine_interface.cpp`` demonstrates the full
  TurbineInterface pipeline, producing NetCDF data that feeds the visualization workflow in
  :doc:`vtk_visualization`.

Integration with VTK Visualization
-----------------------------------

NetCDF state files are the direct input to the ParaView workflow documented in :doc:`vtk_visualization`.
Use ``generate_vtk_output.py`` to convert the desired ``*.nc`` file alongside the matching mesh
connectivity YAML, then open the resulting ``simulation.pvd`` file in ParaView.

Best Practices
--------------

File Management
~~~~~~~~~~~~~~~

- **Unique output paths**: Use descriptive output paths to avoid overwriting results from
  different simulations.
- **File closure**: Call ``Close()`` periodically in long runs (and ``Open()`` before resuming
  writes) to reduce the risk of partial or locked files after unexpected termination.

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

- **Enable buffering**: For simulations with many timesteps (>100), use buffering to reduce
  I/O overhead.
- **Selective components**: Only write state components needed for analysis to reduce file size
  and write time.
- **Selective timesteps**: Gate calls to the writers and emit only every Nth solve if full
  temporal resolution is unnecessary.
- **Chunking**: Chunk sizes are derived from the buffer setting and optimized for time‑sequential
  writes; most workflows do not require manual tuning.
