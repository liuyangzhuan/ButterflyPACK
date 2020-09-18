cd ..
sed -i "" 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/cmake-3.18.2-rik7fyr2qswxr2l7jlyw3zbx2dztrbic/bin/cmake .. \
	-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -finit-real=nan" \
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-DTPL_BLAS_LIBRARIES="/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openblas-0.3.10-gkdwhwee5ga6m26iyurnjlyuz3cqznuq/lib/libopenblas.dylib" \
	-DTPL_LAPACK_LIBRARIES="/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openblas-0.3.10-gkdwhwee5ga6m26iyurnjlyuz3cqznuq/lib/libopenblas.dylib" \
	-DTPL_SCALAPACK_LIBRARIES="/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/netlib-scalapack-2.1.0-erlvlp4bapnc5k2sxbu36entw4vpkl26/lib/libscalapack.dylib" \
	-DTPL_ARPACK_LIBRARIES="/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/arpack-ng-3.7.0-p5kn37ic2b7nz3xbt6fmipz7nvtusqbg/lib/libparpack.dylib;/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/arpack-ng-3.7.0-p5kn37ic2b7nz3xbt6fmipz7nvtusqbg/lib/libarpack.dylib" \
	-DCMAKE_Fortran_COMPILER=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openmpi-3.1.6-c7fikfgh76qcrpjqw7uv4vxsbe5o3yqc//bin/mpif90 \
	-DMPI_Fortran_COMPILER=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openmpi-3.1.6-c7fikfgh76qcrpjqw7uv4vxsbe5o3yqc//bin/mpif90 \
	-DCMAKE_CXX_COMPILER=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openmpi-3.1.6-c7fikfgh76qcrpjqw7uv4vxsbe5o3yqc//bin/mpic++ \
	-DMPI_CXX_COMPILER=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openmpi-3.1.6-c7fikfgh76qcrpjqw7uv4vxsbe5o3yqc//bin/mpic++ \
	-DCMAKE_C_COMPILER=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openmpi-3.1.6-c7fikfgh76qcrpjqw7uv4vxsbe5o3yqc//bin/mpicc \
	-DMPI_C_COMPILER=/Users/liuyangzhuan/Desktop/research/codes/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.0/openmpi-3.1.6-c7fikfgh76qcrpjqw7uv4vxsbe5o3yqc//bin/mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON






	#-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -finit-real=nan" \
	# -DTPL_BLAS_LIBRARIES="" \
	# -DTPL_LAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so" \
	# -DTPL_SCALAPACK_LIBRARIES="/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so" \
