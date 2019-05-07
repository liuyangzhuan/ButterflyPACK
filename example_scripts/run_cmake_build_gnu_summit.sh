module load netlib-lapack/3.8.0
module load gcc/7.4.0
module load cmake/3.11.3
 
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
	-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -finit-real=nan"\
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=OFF \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpiCC\
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/netlib-lapack-3.8.0-wcabdyqhdi5rooxbkqa6x5d7hxyxwdkm/lib64/libblas.so" \
	-DTPL_LAPACK_LIBRARIES="/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/netlib-lapack-3.8.0-wcabdyqhdi5rooxbkqa6x5d7hxyxwdkm/lib64/liblapack.so" \
	-DTPL_SCALAPACK_LIBRARIES="/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/netlib-scalapack-2.0.2-gc7bg7j2zg7hd26sz4ckwcfnwysuqvcd/lib/libscalapack.so"
	# -DTPL_ARPACK_LIBRARIES="/ccs/home/liuyangz/my_software/arpack-ng-gnu/build/lib/libarpack.so;/ccs/home/liuyangz/my_software/arpack-ng-gnu/build/lib/libparpack.so"
# -DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check " \
