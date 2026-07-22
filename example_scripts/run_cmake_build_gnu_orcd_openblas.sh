#!/bin/bash
# Build ButterflyPACK on the MIT ORCD cluster (GNU + OpenMPI + OpenBLAS).
# Adapted from run_cmake_build_gnu_ubuntu.sh.

# ---- modules (hierarchical: load gcc first, then MPI, then libraries) ----
module purge
module load gcc/12.2.0
module load openmpi/4.1.4                 # or openmpi/5.0.8
module load openblas/0.3.26               # provides BLAS + LAPACK
module load netlib-scalapack/2.2.0        # ScaLAPACK (depends on MPI + BLAS/LAPACK)
module load cmake/3.27.9

# ---- library locations ----
# Fill in the install prefixes below. Find each with (after loading the modules):
#   module show openblas/0.3.26          -> look for the ".../openblas-0.3.26-<hash>" prefix
#   module show netlib-scalapack/2.2.0   -> look for the ".../netlib-scalapack-2.2.0-<hash>" prefix
# OpenBLAS ships LAPACK inside libopenblas.so, so BLAS and LAPACK use the same file.
OPENBLAS_LIB="/orcd/software/core/001/spack/pkg/openblas/0.3.26/ro5tivv/lib/libopenblas.so"
SCALAPACK_LIB="/orcd/software/core/001/spack/pkg/netlib-scalapack/2.2.0/ziv7g2h/lib/libscalapack.so"

cd ..
sed -i 's/\r$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt DartConfiguration.tcl CTestTestfile.cmake cmake_install.cmake CMakeFiles
rm -rf SRC_DOUBLE SRC_DOUBLECOMPLEX SRC_SINGLE SRC_COMPLEX

cmake .. \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE -fallow-argument-mismatch" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-Denable_python=OFF \
	-Denable_doc=OFF \
	-DTPL_BLAS_LIBRARIES="${OPENBLAS_LIB}" \
	-DTPL_LAPACK_LIBRARIES="${OPENBLAS_LIB}" \
	-DTPL_SCALAPACK_LIBRARIES="${SCALAPACK_LIB}" \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DOpenMP_C_FLAGS="-fopenmp" \
	-DOpenMP_C_LIB_NAMES="gomp" \
	-DOpenMP_CXX_FLAGS="-fopenmp" \
	-DOpenMP_CXX_LIB_NAMES="gomp" \
	-DOpenMP_Fortran_FLAGS="-fopenmp" \
	-DOpenMP_Fortran_LIB_NAMES="gomp" \
	-DOpenMP_omp_LIBRARY=$(gcc --print-file-name=libgomp.so)

# cvie3d pulls in its dependencies automatically
make cvie3d -j16
