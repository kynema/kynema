# Run this from the root of the Kynema repository e.g.
# docker build -f recipes/Dockerfile -t kynema-fmb:latest .

FROM ubuntu:latest

# Use bash as the default shell
SHELL ["/bin/bash", "-c"]

# Step 1: Install required packages and tools
RUN apt-get update -qq && \
    apt-get install -y \
        build-essential \
        clang-tidy \
        cmake \
        curl \
        g++ \
        gfortran \
        git \
        libcurl4-openssl-dev \
        make \
        python3 \
        python3-pip \
        unzip \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Step 2: Get spack and install necessary libraries using it
RUN git clone https://github.com/spack/spack.git && \
    source spack/share/spack/setup-env.sh && \
    spack compiler find && \
    spack external find && \
    spack mirror add spack-public https://mirror.spack.io && \
    spack install googletest && \
    spack --insecure install cppcheck && \
    spack install yaml-cpp && \
    spack install kokkos && \
    spack install kokkos-kernels && \
    spack install netcdf-c && \
    spack install suite-sparse && \
    spack install eigen && \
    spack install rosco && \
    spack install openfast

# Step 3: Set up the working directory and copy the Kynema repository
WORKDIR /projects/kynema-fmb
COPY . .

# Step 4: Load necessary libraries and build the project
ARG CMAKE_OPTIONS="-DKYNEMA_FMB_WRITE_OUTPUTS=ON \
    -DKYNEMA_FMB_ENABLE_OPENFAST_ADI=ON \
    -DKYNEMA_FMB_ENABLE_ROSCO_CONTROLLER=ON \
    -DCMAKE_BUILD_TYPE=Release"

ARG N_PARALLEL_JOBS=4

RUN . /spack/share/spack/setup-env.sh && \
    spack load kokkos kokkos-kernels netdcf-c suite-sparse eigen googletest yaml-cpp rosco openfast && \
    mkdir docker-build && \
    cd docker-build && \
    cmake .. ${CMAKE_OPTIONS} && \
    cmake --build . -j ${N_PARALLEL_JOBS} && \
    ctest --output-on-failure -j ${N_PARALLEL_JOBS}
