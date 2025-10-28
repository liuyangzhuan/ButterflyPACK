module load PrgEnv-intel
MKL_DIR=$MKLROOT
module load PrgEnv-gnu
module load cmake
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
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="$MKL_DIR//lib/intel64/libmkl_gf_lp64.so;$MKL_DIR//lib/intel64/libmkl_gnu_thread.so;$MKL_DIR//lib/intel64/libmkl_core.so" \
	-DTPL_LAPACK_LIBRARIES="$MKL_DIR//lib/intel64/libmkl_gf_lp64.so;$MKL_DIR//lib/intel64/libmkl_gnu_thread.so;$MKL_DIR//lib/intel64/libmkl_core.so" \
	-DTPL_SCALAPACK_LIBRARIES="$MKL_DIR//lib/intel64/libmkl_blacs_intelmpi_lp64.so;$MKL_DIR//lib/intel64/libmkl_scalapack_lp64.so"

make ctest -j
# make ie2d -j16
# make ie3dport -j16
make ie3d -j16
# make frankben -j16   
# make frankben_t -j16   
# make install

	# -DTPL_ARPACK_LIBRARIES="/global/homes/l/liuyangz/Perlmutter/my_research/arpack-ng-gnu_libsci/build/lib/libparpack.so"