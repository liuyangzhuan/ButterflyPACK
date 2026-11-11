#!/bin/bash
# Build ButterflyPACK and the acoustic Chebyshev driver ie3dacousticcheb (EXAMPLE/Acoustic_SURF_*_CHEB.f90) on Perlmutter CPU
# nodes, with the modules and libraries of run_cmake_build_gnu_perlmutter_openblas.sh. Run it from example_scripts:
#   bash run_cmake_build_gnu_perlmutter_acoustic_cheb.sh           configure from scratch in ../$BUILD_DIR and build
#   bash run_cmake_build_gnu_perlmutter_acoustic_cheb.sh rebuild   after editing SRC/*.f90 or EXAMPLE/*.f90: regenerate the
#                                                                  precision copies (SRC_DOUBLECOMPLEX etc.) and rebuild
# BUILD_DIR (default build) is the build directory under the repository root. cmake runs PrecisionPreprocessing.sh only when
# it configures, so plain "make" after editing SRC/*.f90 compiles the old copies; use the rebuild mode instead.
# The library links the sequential OpenBLAS; runit_acoustic_cheb_perlmutter.sh preloads the OpenMP OpenBLAS and tcmalloc.
module load PrgEnv-gnu
module load cmake
module load python/3.13
module unload cray-libsci
module load cray-fftw

BUILD_DIR=${BUILD_DIR:-build}
cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
export CRAYPE_LINK_TYPE=dynamic
export ZFP_INSTALL_DIR=$CFS/m2957/liuyangz/my_research/zfp-1.0.0_gcc_perlmutter/install
export PARSEC_INSTALL_DIR=/global/cfs/cdirs/m2957/liuyangz/my_software/parsec_pr759/install
export TCMALLOC_ROOT=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/gperftools-2.18.1
export TCMALLOC_LIBRARY="$TCMALLOC_ROOT/lib/libtcmalloc.so"
export TCMALLOC_LINK_FLAGS="-Wl,--no-as-needed $TCMALLOC_LIBRARY -Wl,--as-needed -Wl,-rpath,$TCMALLOC_ROOT/lib"
export CMAKE_PREFIX_PATH="$TCMALLOC_ROOT:$PARSEC_INSTALL_DIR:${CMAKE_PREFIX_PATH:-}"
export PATH="$PARSEC_INSTALL_DIR/bin:${PATH}"
export LD_LIBRARY_PATH="$TCMALLOC_ROOT/lib:$PARSEC_INSTALL_DIR/lib64:${LD_LIBRARY_PATH:-}"

if [ "$1" == "rebuild" ]; then
	bash PrecisionPreprocessing.sh OFF > /dev/null
	cd $BUILD_DIR
	make ie3dacousticcheb -j32
	exit $?
fi

mkdir -p $BUILD_DIR
cd $BUILD_DIR
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
rm -rf SRC_DOUBLE SRC_DOUBLECOMPLEX SRC_SINGLE SRC_COMPLEX
cmake .. \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE" \
	-DCMAKE_CXX_FLAGS="" \
	-DCMAKE_EXE_LINKER_FLAGS="$TCMALLOC_LINK_FLAGS" \
	-DCMAKE_SHARED_LINKER_FLAGS="$TCMALLOC_LINK_FLAGS" \
	-DBUILD_SHARED_LIBS=ON \
	-Denable_python=OFF \
	-Denable_parsec=ON \
	-Denable_fftw=ON \
	-Denable_toplevel_openmp=OFF \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DOpenMP_C_FLAGS="-fopenmp" \
	-DOpenMP_C_LIB_NAMES="gomp" \
	-DOpenMP_Fortran_FLAGS="-fopenmp" \
	-DOpenMP_Fortran_LIB_NAMES="gomp" \
	-DOpenMP_omp_LIBRARY=$(gcc --print-file-name=libgomp.so) \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib64/libzFORp.so;$ZFP_INSTALL_DIR/lib64/libzfp.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
	-DTPL_BLAS_LIBRARIES="$CFS/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib//libopenblas.so" \
	-DTPL_LAPACK_LIBRARIES="$CFS/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib//libopenblas.so" \
	-DTPL_SCALAPACK_LIBRARIES="$CFS/m2957/lib/lib/PrgEnv-gnu/scalapack-2.2.0_sequential/build/install/lib/libscalapack.so"
make ie3dacousticcheb -j32
