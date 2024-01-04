module load PrgEnv-gnu
module load cmake

cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
    -DCMAKE_INSTALL_LIBDIR=./lib \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE" \
    -DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu_82_mpi_mp.so" \
	-DTPL_LAPACK_LIBRARIES="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu_82_mpi_mp.so" \
	-DTPL_SCALAPACK_LIBRARIES="${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu_82_mpi_mp.so" \
	-DTPL_ARPACK_LIBRARIES="/global/homes/l/liuyangz/Perlmutter/my_research/arpack-ng-gnu_libsci/build/lib/libparpack.so"
make ie2d -j16
make ie3dport -j16
make frankben -j16   
make frankben_t -j16   
make install
