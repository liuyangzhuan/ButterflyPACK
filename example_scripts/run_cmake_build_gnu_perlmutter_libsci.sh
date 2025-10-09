module load PrgEnv-gnu
module load cmake

LIBSCI_DIR=$CRAY_LIBSCI_PREFIX
module unload cray-libsci

cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
export ZFP_INSTALL_DIR=$CFS/m2957/liuyangz/my_research/zfp-1.0.0_gcc_perlmutter/install

rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
    -DCMAKE_INSTALL_LIBDIR=./lib \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE" \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib64/libzFORp.so;$ZFP_INSTALL_DIR/lib64/libzfp.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
    -DCMAKE_CXX_FLAGS="" \
	-Denable_toplevel_openmp=OFF \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DOpenMP_C_FLAGS="-fopenmp" \
	-DOpenMP_C_LIB_NAMES="gomp" \
	-DOpenMP_Fortran_FLAGS="-fopenmp" \
	-DOpenMP_Fortran_LIB_NAMES="gomp" \
	-DOpenMP_omp_LIBRARY=$(gcc --print-file-name=libgomp.so) \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="$LIBSCI_DIR/lib/libsci_gnu_mpi_mp.so" \
	-DTPL_LAPACK_LIBRARIES="$LIBSCI_DIR/lib/libsci_gnu_mpi_mp.so" \
	-DTPL_SCALAPACK_LIBRARIES="$LIBSCI_DIR/lib/libsci_gnu_mpi_mp.so"

make ctest -j
# make ie2d -j16
# make ie3dport -j16
make ie3d -j16
# make frankben -j16   
# make frankben_t -j16   
# make install

	# -DTPL_ARPACK_LIBRARIES="/global/homes/l/liuyangz/Perlmutter/my_research/arpack-ng-gnu_libsci/build/lib/libparpack.so"