#!/bin/bash
# Perlmutter GPU build: sequential OpenBLAS (no LibSci), MPI, OpenMP, PaRSEC,
# and the GPU backend of the H2 Color factorization (enable_h2_gpu, MAGMA on
# A100).  Configures into ../build_gpu, next to the CPU build in ../build.
module load PrgEnv-gnu cray-fftw cmake python cudatoolkit craype-accel-nvidia80
module unload cray-libsci

cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build_gpu
cd build_gpu
export CRAYPE_LINK_TYPE=dynamic

LIBDIR=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu
MAGMA_DIR=/global/cfs/cdirs/m2957/lib/magma_master
ZFP_INSTALL_DIR=/global/cfs/cdirs/m2957/liuyangz/my_research/zfp-1.0.0_gcc_perlmutter/install

rm -rf CMakeCache.txt CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-Denable_mpi=ON \
	-Denable_openmp=ON \
	-Denable_toplevel_openmp=OFF \
	-Denable_fftw=ON \
	-Denable_parsec=ON \
	-Denable_python=ON \
	-Denable_h2_gpu=ON \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DPaRSEC_DIR=/global/cfs/cdirs/m2957/liuyangz/my_software/parsec_pr759/install/share/cmake/parsec \
	-DTPL_BLAS_LIBRARIES="$LIBDIR/OpenBLAS_sequential/build/install/lib/libopenblas.so" \
	-DTPL_LAPACK_LIBRARIES="$LIBDIR/OpenBLAS_sequential/build/install/lib/libopenblas.so" \
	-DTPL_SCALAPACK_LIBRARIES="$LIBDIR/scalapack-2.2.0_sequential/build/install/lib/libscalapack.so" \
	-DTPL_FFTW_LIBRARIES="/opt/cray/pe/fftw/3.3.10.11/x86_milan/lib/libfftw3.so;/opt/cray/pe/fftw/3.3.10.11/x86_milan/lib/libfftw3f.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib64/libzFORp.so;$ZFP_INSTALL_DIR/lib64/libzfp.so" \
	-DTPL_H2_MAGMA_INCLUDE_DIRS="$MAGMA_DIR/include" \
	-DTPL_H2_MAGMA_LIBRARIES="$MAGMA_DIR/lib/libmagma.so"

make -j 32 claplace3d_h2 cmatern1d_h2 cmatern2d_h2 cvie3d_h2
