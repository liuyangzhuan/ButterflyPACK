spackdir=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-bigsur-x86_64/apple-clang-14.0.0
cd ..
sed -i "" 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
${spackdir}/cmake-3.21.2-4m4wepkygdb3trocubzle24hsrh4mddv/bin/cmake .. \
	-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -finit-real=nan" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-DTPL_BLAS_LIBRARIES="${spackdir}/openblas-0.3.17-cs3jktygp27xfmvkdv5uh6iw5tvhy2xz/lib/libopenblas.dylib" \
	-DTPL_LAPACK_LIBRARIES="${spackdir}/openblas-0.3.17-cs3jktygp27xfmvkdv5uh6iw5tvhy2xz/lib/libopenblas.dylib" \
	-DTPL_SCALAPACK_LIBRARIES="${spackdir}/netlib-scalapack-2.1.0-eztdodlcfmw4qugmbjr2eybes2q6om3q/lib/libscalapack.dylib" \
	-DTPL_ARPACK_LIBRARIES="${spackdir}/arpack-ng-3.8.0-qtmg4htq24esntjtlt5lccwipdsheexn/lib/libparpack.dylib;${spackdir}/arpack-ng-3.8.0-qtmg4htq24esntjtlt5lccwipdsheexn/lib/libarpack.dylib" \
	-DCMAKE_Fortran_COMPILER=${spackdir}/openmpi-4.1.1-mmdjfqrs5cn73jw5dzmfmo6nknspctt3/bin/mpif90 \
	-DMPI_Fortran_COMPILER=${spackdir}/openmpi-4.1.1-mmdjfqrs5cn73jw5dzmfmo6nknspctt3/bin/mpif90 \
	-DCMAKE_CXX_COMPILER=${spackdir}/openmpi-4.1.1-mmdjfqrs5cn73jw5dzmfmo6nknspctt3/bin/mpic++ \
	-DMPI_CXX_COMPILER=${spackdir}/openmpi-4.1.1-mmdjfqrs5cn73jw5dzmfmo6nknspctt3/bin/mpic++ \
	-DCMAKE_C_COMPILER=${spackdir}/openmpi-4.1.1-mmdjfqrs5cn73jw5dzmfmo6nknspctt3/bin/mpicc \
	-DMPI_C_COMPILER=${spackdir}/openmpi-4.1.1-mmdjfqrs5cn73jw5dzmfmo6nknspctt3/bin/mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON






	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
