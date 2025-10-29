# -----------------------------------------------------------------------------
# Run this from the root of the Kynema repository with the following arguments:
# ./recipes/build_recipe_macos.sh <path_to_kynema_root> <path_to_spack_root>
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Note that many systems have ROCm already built and properly configured.  If this
# is the case for you, make sure that all modules are loaded and paths set correctly.
# You may have to manually set some of this for Spack configuration.
# You may also need to add the amdgpu_target option for your platform to the Spack
# spec.
# -----------------------------------------------------------------------------

# Exit on error
set -e

# Check if FC (Fortran Compiler) is set
if [ -z "${FC}" ]; then
    echo "Error: FC (Fortran Compiler) environment variable is not set."
    echo "Please set FC to your preferred Fortran compiler before running this script."
    echo "For example: export FC=/path/to/your/gfortran"
    exit 1
fi

# Clone Kynema repository if not present in the provided path
kynema_path=$(dirname "$1") # Get the directory of the provided path
if [ ! -d "$kynema_path/kynema" ]; then
    git clone --recursive https://github.com/kynema/kynema.git $kynema_path/kynema
fi
cd $kynema_path/kynema

# Install Spack if not present in the provided path
spack_path=$(dirname "$2") # Get the directory of the provided path
if [ ! -d "$spack_path/spack" ]; then
    git clone https://github.com/spack/spack.git $spack_path/spack
fi
source $spack_path/spack/share/spack/setup-env.sh
spack compiler find
spack external find

# Install required packages using Spack if not present
install_if_missing() {
    if ! spack find -l "$1" | grep -q "$1"; then
        echo "Installing $1..."
        spack install "$1"
    else
        echo "$1 is already installed."
    fi
}

install_if_missing googletest
install_if_missing yaml-cpp
install_if_missing kokkos+rocm
install_if_missing kokkos-kernels+rocsparse+rocblas
install_if_missing netcdf-c
install_if_missing suite-sparse
install_if_missing eigen
install_if_missing rosco
install_if_missing openfast
#install_if_missing cppcheck # add if CppCheck is needed
#install_if_missing llvm # add if clang-tidy is needed

spack load kokkos kokkos-kernels netdcf-c suite-sparse eigen googletest yaml-cpp rosco openfast #llvm cppcheck

# Build Kynema with the specified options
mkdir -p build-from-script
cd build-from-script
cmake .. \
  -DCMAKE_CXX_COMPILER=hipcc \
  -DKynema_WRITE_OUTPUTS=ON \
  -DKynema_ENABLE_OPENFAST_ADI=ON \
  -DKynema_ENABLE_ROSCO_CONTROLLER=ON \
  -DKynema_ENABLE_KLU=ON \
  -DCMAKE_BUILD_TYPE="Release"

cmake --build .
