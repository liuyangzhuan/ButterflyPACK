export CRAYPE_LINK_TYPE=dynamic
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-DPRNTlevel=2 -DGNU -O0 -pg -cpp -fbacktrace -fimplicit-none -fbounds-check -ffree-line-length-none  -ffixed-line-length-none -fopenmp -Wconversion -lpthread " \
	-DCMAKE_CXX_FLAGS="-std=c++11 -O3 -fopenmp -lifcore" \
	-DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	-DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-Denable_complex=OFF
	# -DCMAKE_Fortran_FLAGS="-cpp -DGNU -DPRNTlevel=0 -O3 -ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -ffree-line-length-none -ffixed-line-length-none -fopenmp " \
	# -DCMAKE_Fortran_FLAGS="-DPRNTlevel=2 -DGNU -O0 -pg -cpp -fbacktrace -ffpe-trap=zero,overflow,underflow -fimplicit-none -fbounds-check -ffree-line-length-none  -ffixed-line-length-none -fopenmp -Wconversion -lpthread " \	
	
