module purge
# module load pgi
# module load scalapack-netlib/pgi-19.10/2.0.2
module load PrgEnv-pgi

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
	-DCMAKE_Fortran_FLAGS="" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=OFF \
	-DTPL_BLAS_LIBRARIES="/opt/pgi/linux86-64-llvm/19.10/lib/libomp.so;/opt/pgi/linux86-64-llvm/19.10/lib/libblas.so" \
	-DTPL_LAPACK_LIBRARIES="/opt/pgi/linux86-64-llvm/19.10/lib/liblapack.so" \
	-DTPL_SCALAPACK_LIBRARIES="/opt/pgi/linux86-64-llvm/19.10/lib/scalapack/scalapack-2.0.2/openmpi-3.1.3/lib/libscalapack.a" \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON






	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
