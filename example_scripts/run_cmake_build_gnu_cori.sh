module swap PrgEnv-intel PrgEnv-gnu
cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
bash PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
export MKLROOT=/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS=""\
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=OFF \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="" \
	-DTPL_LAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_gf_lp64.so;${MKLROOT}/lib/intel64/libmkl_gnu_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so" \
	-DTPL_SCALAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so;${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so" \
	-DTPL_ARPACK_LIBRARIES="/global/homes/l/liuyangz/Cori/my_software/arpack-ng-gnu/build/lib/libarpack.so;/global/homes/l/liuyangz/Cori/my_software/arpack-ng-gnu/build/lib/libparpack.so"
# -DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check " \
