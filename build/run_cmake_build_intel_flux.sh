cd ..
sh PrecisionPreprocessing.sh
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-cpp -DIntel -nologo -fpe0 -traceback -debug full" \
	-DCMAKE_CXX_FLAGS="-std=c++11 -O3 -qopenmp -qopt-matmul -lifcore " \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=mpifort \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_LAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_intel_lp64.so;${MKLROOT}/lib/intel64/libmkl_intel_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so;-liomp5 -lpthread -lm -ldl" \
	-DTPL_SCALAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so;${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so" 	
	#-DCMAKE_Fortran_FLAGS="-O3 -no-prec-div -axAVX,SSE4.2 -msse2 -align records -parallel -lpthread -I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include/intel64/lp64" \
	#-DCMAKE_CXX_FLAGS="-std=c++11 -O3 -qopt-matmul" \
	# -DCMAKE_Fortran_FLAGS="-nologo -fpe0 -traceback -debug full -debug parallel -O0 -g -check bounds -qopenmp -parallel -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64" \
	# -DCMAKE_CXX_FLAGS="-O0 -g -std=c++11 -qopenmp -debug parallel -traceback" \
	# -DTPL_LAPACK95_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blas95_lp64.a;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a" \
	# -DTPL_LAPACK95_INCLUDE_DIRS="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include/intel64/lp64" \
#	-DCMAKE_Fortran_FLAGS="-cpp -DIntel -no-prec-div -align records -parallel -lpthread" \
#	-DCMAKE_Fortran_FLAGS="-cpp -DIntel -nologo -fpe0 -traceback -debug full" \
