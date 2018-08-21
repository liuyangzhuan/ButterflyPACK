export CRAYPE_LINK_TYPE=dynamic
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_CXX_FLAGS="-std=c++11" \
	-DCMAKE_Fortran_FLAGS="-cpp -DIntel -DPRNTlevel=0 -O3 -no-prec-div -axAVX,SSE4.2 -msse2 -align records -parallel -lpthread" \
	-DTPL_SCALAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so;${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so" \
	-DBUILD_SHARED_LIBS=OFF \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
        -DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-Denable_complex=OFF
	# -DCMAKE_Fortran_FLAGS="-cpp -DIntel -DPRNTlevel=0 -nologo -traceback -fpe0 -debug full -O0 -check bounds" \ 
	# -DCMAKE_Fortran_FLAGS="-cpp -DIntel -DPRNTlevel=0 -O3 -no-prec-div -axAVX,SSE4.2 -msse2 -align records -parallel -lpthread" \
	
	
	