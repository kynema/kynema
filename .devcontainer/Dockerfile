
FROM debian:testing

ENV TZ=America/New_York

# Development tools/libraries
RUN apt update -y && \
    apt install -y git git-lfs wget nodejs \
    gcc g++ gfortran cmake gdb valgrind \
    libopenblas-serial-dev libopenmpi-dev \
    libnetcdf-dev libhdf5-dev libhdf5-serial-dev
