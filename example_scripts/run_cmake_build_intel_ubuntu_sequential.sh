module purge
module load compiler
module load mpi
module load mkl
module load cmake

export ZFP_INSTALL_DIR=/home/administrator/Desktop/Research/zfp_intel/install/


cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
# bash PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="" \
	-Denable_mpi=OFF \
	-DCMAKE_CXX_FLAGS="-I$MKLROOT/include " \
	-DTPL_BLAS_LIBRARIES="$MKLROOT//lib/intel64/libmkl_intel_lp64.so;$MKLROOT//lib/intel64/libmkl_intel_thread.so;$MKLROOT//lib/intel64/libmkl_core.so" \
	-DTPL_LAPACK_LIBRARIES="$MKLROOT//lib/intel64/libmkl_intel_lp64.so;$MKLROOT//lib/intel64/libmkl_intel_thread.so;$MKLROOT//lib/intel64/libmkl_core.so" \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib/libzFORp.so;$ZFP_INSTALL_DIR/lib/libzfp.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
	-DBUILD_SHARED_LIBS=OFF \
	-DCMAKE_Fortran_COMPILER=ifort \
	-DCMAKE_CXX_COMPILER=icpc \
	-DCMAKE_C_COMPILER=icc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
	# -DCMAKE_Fortran_FLAGS="-no-prec-div -align records -parallel -lpthread" \
make ie2d
make ctest_simple_seq_newapi
	# -DCMAKE_EXE_LINKER_FLAGS="-qopenmp -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core" \
