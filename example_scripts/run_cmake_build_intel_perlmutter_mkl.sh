module load PrgEnv-intel
module load cmake
module unload cray-libsci



cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
# bash PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
export ZFP_INSTALL_DIR=$CFS/m2957/liuyangz/my_research/zfp-1.0.0_intel_perlmutter/install


rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="" \
	-DCMAKE_CXX_FLAGS="-I$MKLROOT/include " \
	-DTPL_BLAS_LIBRARIES="$MKLROOT//lib/intel64/libmkl_intel_lp64.so;$MKLROOT//lib/intel64/libmkl_intel_thread.so;$MKLROOT//lib/intel64/libmkl_core.so" \
	-DTPL_LAPACK_LIBRARIES="$MKLROOT//lib/intel64/libmkl_intel_lp64.so;$MKLROOT//lib/intel64/libmkl_intel_thread.so;$MKLROOT//lib/intel64/libmkl_core.so" \
	-DTPL_SCALAPACK_LIBRARIES="$MKLROOT//lib/intel64/libmkl_blacs_intelmpi_lp64.so;$MKLROOT//lib/intel64/libmkl_scalapack_lp64.so" \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib/libzFORp.so;$ZFP_INSTALL_DIR/lib/libzfp.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DMPI_CXX_COMPILER=CC \
	-DMPI_C_COMPILER=cc \
	-DMPI_Fortran_COMPILER=ftn \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
	# -DCMAKE_Fortran_FLAGS="-no-prec-div -align records -parallel -lpthread" \
make ie3d -j
	# -DCMAKE_EXE_LINKER_FLAGS="-qopenmp -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core" \
