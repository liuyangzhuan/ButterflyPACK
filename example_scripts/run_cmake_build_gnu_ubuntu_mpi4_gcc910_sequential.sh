module purge
module load gcc/9.1.0
# module load openmpi/gcc-9.1.0/4.0.1
# module load scalapack-netlib/gcc-9.1.0/2.2.0
module load cmake/3.19.2
module unload python
export PATH=/home/administrator/Desktop/Research/GPTune_master/env/bin/:$PATH

export ZFP_INSTALL_DIR=/home/administrator/Desktop/Research/zfp/install/

cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
# bash PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -finit-real=nan" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-Denable_doc=OFF \
	-Denable_openmp=ON \
	-Denable_mpi=OFF \
	-Denable_python=ON \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib/libzFORp.so;$ZFP_INSTALL_DIR/lib/libzfp.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
	-DTPL_BLAS_LIBRARIES="/usr/lib/x86_64-linux-gnu/libblas.so" \
	-DTPL_LAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/liblapack.so" \
	-DCMAKE_Fortran_COMPILER=$F90 \
	-DCMAKE_CXX_COMPILER=$CXX \
	-DCMAKE_C_COMPILER=$CC \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_INSTALL_LIBDIR=./lib \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON

make krr_seq -j16
make ctest_simple_seq -j16
make butterflypack_cpp
make ctest_simple_seq_newapi
make install -j


	# -DTPL_ARPACK_LIBRARIES="/home/administrator/Desktop/Software/arpack-ng/build/lib/libarpack.so;/home/administrator/Desktop/Software/arpack-ng/build/lib/libparpack.so" \


#	-DTPL_SCALAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/libscalapack.so" \


	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
	# -DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib/libzFORp.so;$ZFP_INSTALL_DIR/lib/libzfp.so" \
	# -DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
