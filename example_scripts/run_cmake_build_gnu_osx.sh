cd ..
sed -i "" 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/cmake-3.15.1-fg3vqdb6wqkb7vwi7vdlnvi5qhuyvd34/bin/cmake .. \
	-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -finit-real=nan" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=OFF \
	-DTPL_BLAS_LIBRARIES="/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openblas-0.3.7-pmh2e7ebbmrdk65p34rpqyb3lhlbduo5/lib/libopenblas.dylib" \
	-DTPL_LAPACK_LIBRARIES="" \
	-DTPL_SCALAPACK_LIBRARIES="/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/netlib-scalapack-2.0.2-3xcgrjzdj32awb5fifhed3bn6n7qd3hh/lib/libscalapack.dylib" \
	-DTPL_ARPACK_LIBRARIES="/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/arpack-ng-3.7.0-tp6enbmmgupgdzyacydjbfvykqhynmnp/lib/libparpack.dylib;/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/arpack-ng-3.7.0-tp6enbmmgupgdzyacydjbfvykqhynmnp/lib/libarpack.dylib" \
	-DCMAKE_Fortran_COMPILER=/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openmpi-3.1.4-yyoaovisks4oywyfecxlz2gl4v5mrnp2/bin/mpif90 \
	-DMPI_Fortran_COMPILER=/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openmpi-3.1.4-yyoaovisks4oywyfecxlz2gl4v5mrnp2/bin/mpif90 \
	-DCMAKE_CXX_COMPILER=/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openmpi-3.1.4-yyoaovisks4oywyfecxlz2gl4v5mrnp2/bin/mpic++ \
	-DMPI_CXX_COMPILER=/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openmpi-3.1.4-yyoaovisks4oywyfecxlz2gl4v5mrnp2/bin/mpic++ \
	-DCMAKE_C_COMPILER=/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openmpi-3.1.4-yyoaovisks4oywyfecxlz2gl4v5mrnp2/bin/mpicc \
	-DMPI_C_COMPILER=/Users/yangliu/Desktop/my_software/spack/opt/spack/darwin-highsierra-x86_64/clang-10.0.0-apple/openmpi-3.1.4-yyoaovisks4oywyfecxlz2gl4v5mrnp2/bin/mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON






	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
