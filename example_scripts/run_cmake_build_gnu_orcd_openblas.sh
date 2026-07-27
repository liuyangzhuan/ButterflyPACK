#!/bin/bash
# Build ButterflyPACK on the MIT ORCD cluster (GNU + OpenMPI + OpenBLAS).
# Adapted from run_cmake_build_gnu_ubuntu.sh.

# ---- modules (hierarchical: load gcc first, then MPI, then libraries) ----
module purge
module load gcc/12.2.0
module load openmpi/4.1.4                 # or openmpi/5.0.8
# module load openblas/0.3.26             # too old (no cblas_?gemm_batch); using $HOME/OpenBLAS build instead
module load netlib-scalapack/2.2.0        # ScaLAPACK (depends on MPI + BLAS/LAPACK)
module load fftw/3.3.10                    # FFTW (serial, double) needed by kernel.hpp
module load cmake/3.27.9

# ---- library locations ----
# OpenBLAS: built from source in $HOME/OpenBLAS (v0.3.29+, includes cblas_?gemm_batch,
# which the module 0.3.26 lacked). OpenBLAS ships LAPACK inside libopenblas.so, so
# BLAS and LAPACK point to the same file.
OPENBLAS_LIB="$HOME/OpenBLAS/libopenblas.so"
SCALAPACK_LIB="/orcd/software/core/001/spack/pkg/netlib-scalapack/2.2.0/ziv7g2h/lib/libscalapack.so"

# FFTW: CMake has no FFTW logic, so we add its include (compile) and lib (link) by hand.
# Find the prefix after `module load fftw/3.3.10` via:  env | grep -i fftw
FFTW_PREFIX="/orcd/software/core/001/spack/pkg/fftw/3.3.10/2qziucy"

# ZFP: built from source into $HOME/zfp-install (git clone LLNL/zfp; cmake -DBUILD_ZFORP=ON;
# make install). Providing TPL_ZFP_* makes CMake define HAVE_ZFP so the *_ZFP_Compress
# imports in CPP_INTERFACE/m_BPACK_utilities.f90 resolve.
ZFP_PREFIX="$HOME/zfp-install"

cd ..
sed -i 's/\r$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt DartConfiguration.tcl CTestTestfile.cmake cmake_install.cmake CMakeFiles
rm -rf SRC_DOUBLE SRC_DOUBLECOMPLEX SRC_SINGLE SRC_COMPLEX

# Build type: RelWithDebInfo compiles -O2 -g (keeps optimization, adds debug symbols
# so gdb/UCX backtraces show file:line). Toggle by swapping which line is commented.
#BUILD_TYPE=Release
BUILD_TYPE=RelWithDebInfo

cmake .. \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE -fallow-argument-mismatch" \
	-DCMAKE_CXX_FLAGS="-I${FFTW_PREFIX}/include" \
	-DCMAKE_EXE_LINKER_FLAGS="-L${FFTW_PREFIX}/lib -lfftw3" \
	-DBUILD_SHARED_LIBS=ON \
	-Denable_python=OFF \
	-Denable_doc=OFF \
	-DTPL_BLAS_LIBRARIES="${OPENBLAS_LIB}" \
	-DTPL_LAPACK_LIBRARIES="${OPENBLAS_LIB}" \
	-DTPL_SCALAPACK_LIBRARIES="${SCALAPACK_LIB}" \
	-DTPL_ZFP_LIBRARIES="${ZFP_PREFIX}/lib64/libzFORp.so;${ZFP_PREFIX}/lib64/libzfp.so" \
	-DTPL_ZFP_INCLUDE="${ZFP_PREFIX}/include" \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
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
