cd ..
sh PrecisionPreprocessing.sh
cd build
export CRAYPE_LINK_TYPE=dynamic
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-cpp -DIntel -DPRNTlevel=2 -nologo -traceback -fpe0 -debug full -O0 -check bounds" \
	-DCMAKE_CXX_FLAGS="-std=c++11 -O3 -qopenmp -qopt-matmul -lifcore -I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include " \
	-DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_lp64.a;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.a;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.a;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.a" \
	-DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.a;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.a" \
	-DCMAKE_EXE_LINKER_FLAGS=" -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core" \
	-DBUILD_SHARED_LIBS=OFF \
	-DCMAKE_Fortran_COMPILER=mpiifort \
	-DCMAKE_CXX_COMPILER=mpiicpc \
	-DCMAKE_C_COMPILER=mpiicc \
	-DMPI_CXX_COMPILER=/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpiicpc \
	-DMPI_C_COMPILER=/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpiicc \
	-DMPI_Fortran_COMPILER=/opt/intel/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpiifort \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
	# -Denable_complex=OFF
	# -DCMAKE_Fortran_FLAGS="-cpp -DIntel -DPRNTlevel=0 -nologo -traceback -fpe0 -debug full -O0 -check bounds" \ 
	# -DCMAKE_Fortran_FLAGS="-cpp -DIntel -DPRNTlevel=0 -O3 -no-prec-div -align records -parallel -lpthread" \
