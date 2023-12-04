#!/bin/zsh

source /usr/local/Cellar/modules/4.3.0/init/zsh

module load gcc/9.2.0
module load openmpi/gcc-9.2.0/4.0.2
module load scalapack/2.1.0


cd ..
sed -i "" 's/^M$//' PrecisionPreprocessing.sh
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
	-DTPL_BLAS_LIBRARIES="/usr/local/Cellar/openblas/0.3.7/lib/libblas.dylib" \
	-DTPL_LAPACK_LIBRARIES="/usr/local/Cellar/openblas/0.3.7/lib/liblapack.dylib" \
	-DTPL_SCALAPACK_LIBRARIES="/usr/local/Cellar/scalapack/2.1.0/lib/libscalapack.dylib" \
	-DCMAKE_Fortran_COMPILER=$MPIF90 \
	-DMPI_Fortran_COMPILER=/$MPIF90 \
	-DCMAKE_CXX_COMPILER=$MPICXX \
	-DMPI_CXX_COMPILER=$MPICXX \
	-DCMAKE_C_COMPILER=$MPICC \
	-DMPI_C_COMPILER=$MPICC \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON






	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
