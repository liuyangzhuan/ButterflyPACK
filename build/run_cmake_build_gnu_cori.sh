module swap PrgEnv-intel PrgEnv-gnu
cd ..
sh PrecisionPreprocessing.sh
cd build
export CRAYPE_LINK_TYPE=dynamic
export MKLROOT=/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-DGNU -O0 -pg -cpp -fbacktrace -fimplicit-none -fbounds-check -ffree-line-length-none  -ffixed-line-length-none -fopenmp -Wconversion -lpthread"\
	-DCMAKE_CXX_FLAGS="-std=c++11 -O3 -fopenmp" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_LAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_gf_lp64.so;${MKLROOT}/lib/intel64/libmkl_gnu_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so" \
	-DTPL_SCALAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so;${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so" 	
	# -DCMAKE_Fortran_FLAGS="-DGNU -O3 -cpp -ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -ffree-line-length-none -ffixed-line-length-none -fopenmp " \
	# -DCMAKE_Fortran_FLAGS="-DGNU -O0 -pg -cpp -fbacktrace -fimplicit-none -fbounds-check -ffree-line-length-none  -ffixed-line-length-none -fopenmp -Wconversion -lpthread " \
