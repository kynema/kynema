Introduction
============

Kynema-FMB is a speed-oriented, performance portable library for solving structural dynamics problems.

Code Structure
--------------

Kynema-FMB provides three element types: high-order beam finite elements, mass elements, and geometrically nonlinear spring elements.
When combined with constraints, these elements can be used to describe any number of structural dynamics configurations.
Kynema-FMB's "low level" API is based around the Model class, where the user will manually specify each node and constraint and element memberships.
Additionally, using the low level API gives the user full flexibility in controlling the time stepping proceedure and copying memory to/from device.
For more details on the low level API, see the Heavy Top, Spring-Mass, and Three Blade example problems.

Kynema-FMB also provides a number of "high level" APIs for defining and driving specific problem configurations that are in high demand by our users.
This approach allows users to avoid much of the tedious book keeping needed to specify and drive these common problems.
For more details on the currently provided low level APIs, see the Rigid Body with Three Springs example and the IEA15MW Turbine example.


Kokkos and Kokkos Kernels
-------------------------

Kynema-FMB is built on top of the Kokkos ecosystem, which allows a single code base to target serial CPU, multi-threaded CPU, and GPU platforms from NVIDIA, AMD, and Intel.

Kokkos provides the basic building blocks for expressing algorithms in a performance portable way.
Kokkos' Views wrap the details of memory allocation across platforms and abstract memory access patterns for optimal caching/coalescing behavior on CPU/GPU platforms. 
Additionally, Kokkos' parallel_for construct abstracts the details of kernel launching on different platforms, allowing a single piece of code to represent a simple for loop (when configured for serial execution), an OpenMP accelerated for loop with thread local variables (when configured for OpenMP execution) or a Cuda/HIP kernel launch with tuned grid sizing and shared memory usaged.
Furthermore, Kokkos allows us to express heirarchical parallelism, which is essential for extracting maximum performance on GPU devices when working with high order elements.
Finally, Kokkos provides a SIMD interface, which allows us to manually specify outer-loop SIMD parallelism in key algorithms.

Kokkos Kernels is a library built on top of Kokkos which provides performance-portable BLAS algorithms. 
Kynema-FMB primarily uses Kokkos Kernels' batched algorithms to perform optimized matrix math in its system assembly algorithms.
Kokkos Kernels also provides a sparse matrix library, which is used for manipulating Kynema-FMB's sparse system matrix.

Linear Solvers
--------------

Much of Kynema-FMB's performance is derived from its use of highly optimized sparse direct solvers.
When solving the system matrix that is assembled as part of our time stepping processs, sparse direct solvers combine the robustness of traditional LU factorization with the efficiency of a sparse matrix representation.
The actual peformance of these algorithms is highly dependent on the fill-in of the matrix's sparsity pattern, the built-in reordering algorithm, and the details of the exact solution proceedure.

On CPU, we provide many choices for solvers, since different solvers have their own performance properties for different sparse matrix structures.
Based on our limited testing, SuiteSparse's KLU solver gives the best performance for our problems of interest, but others are provided in case users find otherwise for their application.
Because of this, KLU has recieved the majority of testing - our other CPU-based solvers should be considered experimental.
Our CPU-based sparse linear solvers are as follows:

- KLU: Provided by SuiteSparse, this solver is optimized for circuit simulations, but also gives supurb performance for our problems.
- UMFPACK: A more general solver from SuiteSparse.
- SuperLU: Another popular sparse direct solver, tends to give the second best performance after KLU.
- SuperLU-MT: A multi-threaded version of SuperLU.
- oneMKL PARDISO: Another multi-threaded sparse direct solver from Intel's oneMKL library.  Note that it currently only works for MKL versions 2022 and earlier.

When running on GPU, if you do not enable a GPU-based sparse linear solver, we will automatically copy the sparse system to host and use a configured CPU-based linear system.
For some problems, this may actually be the most performant option, but we also provide sparse direct solvers that run natively on the GPU for problems where that option is faster.
Which technique to employ is likely problem dependent and a topic of ongoing investigation.

On NVIDIA GPUs, we provide two options for GPU-based solvers as follows:

- cuSolverSP: CUDA's legacy native solver, which is currently deprecated.  It is a part of the cuSolver library, so it should be provided with most CUDA installs.  For small problems, this solver is relatively slow, but it scales well with problem size.
- cuDSS: CUDA's next generation sparse direct solver, still in pre-release at the time of this wrting.  This offers superior performance to cuSolverSP and is often faster than CPU-based alternatives.  It is under ongoing development by NVIDIA, so our usage and recommendations are expected to evolve along side it.

On other GPU platforms (such as AMD's ROCm), we do not offer GPU-native solvers due to a lack of vendor provided libraries.
As these platforms' ecosystems evolve, we will investigate any solvers that become available.
For now, it is recommended to use KLU when running on ROCm devices - the performance, even with the round-trip data transfer should be sufficient.

We have not yet explored other linear solver techniques, such as traditional iterative solvers or algebraic multigrid methods.
The relatively low fill-in of the matrices resulting from our structural dynamics problems, as well as our emperical performance results thus far, implies that sparse direct methods will provide superior performance, but we may investigate this explicitly in the future, especially on GPU.
