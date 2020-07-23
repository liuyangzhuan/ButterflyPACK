module purge

module load gcc/7.4.0
module load openmpi/gcc-7.4.0/4.0.1
module load python/gcc-7.4.0/3.7.1
module load scalapack-netlib/gcc-7.4.0/2.0.2
module load cuda/gcc-7.4.0/10.1
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
	-DCMAKE_C_FLAGS="-I/home/administrator/Desktop/Software/magma-2.5.3/install/include -I/usr/local/cuda-10.1/include" \
	-DBUILD_SHARED_LIBS=ON \
	-DTPL_BLAS_LIBRARIES="/usr/lib/x86_64-linux-gnu/libblas.so" \
	-DTPL_LAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/liblapack.so" \
	-DTPL_SCALAPACK_LIBRARIES="/home/administrator/Desktop/Software/scalapack-2.0.2-gcc-7.4.0/build/lib/libscalapack.so" \
	-DTPL_MAGMA_LIBRARIES="/home/administrator/Desktop/Software/magma-2.5.3/install/lib/libmagma.so" \
	-DCMAKE_Fortran_COMPILER=$MPIF90 \
	-DCMAKE_CXX_COMPILER=$MPICXX \
	-DCMAKE_C_COMPILER=$MPICC \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=RelWithDebInfo \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON






	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
