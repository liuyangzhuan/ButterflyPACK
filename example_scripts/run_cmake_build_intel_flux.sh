cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
bash PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="" \
	-DCMAKE_CXX_FLAGS="-lifcore" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=mpifort \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="" \
	-DTPL_LAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_intel_lp64.so;${MKLROOT}/lib/intel64/libmkl_intel_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so;-liomp5 -lpthread -lm -ldl" \
	-DTPL_SCALAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so;${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so"
#	-DCMAKE_Fortran_FLAGS="-no-prec-div -align records -parallel -lpthread" \

