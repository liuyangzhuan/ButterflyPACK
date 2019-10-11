module swap xl gcc/9.1.0
module load mercurial
module load cmake
module load essl
module load cuda
module load netlib-lapack
module load netlib-scalapack
 
cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
# bash PrecisionPreprocessing.sh
mkdir -p build
cd build
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
	-DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -finit-real=nan -ffree-line-length-none -ffixed-line-length-none"\
	-DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=OFF \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpiCC\
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=RelWithDebInfo \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
        -DTPL_BLAS_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so;${OLCF_NETLIB_LAPACK_ROOT}/lib64/libblas.so" \
        -DTPL_LAPACK_LIBRARIES="${OLCF_ESSL_ROOT}/lib64/libessl.so;${OLCF_NETLIB_LAPACK_ROOT}/lib64/liblapack.so" \
        -DTPL_SCALAPACK_LIBRARIES="${OLCF_NETLIB_SCALAPACK_ROOT}/lib/libscalapack.so" \
	-DMPI_CXX_COMPILER=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-9.1.0/spectrum-mpi-10.3.0.1-20190611-2juhkwlddpdydl2yasmuia2slpev6fl5/bin/mpiCC \
	-DMPI_Fortran_COMPILER=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-9.1.0/spectrum-mpi-10.3.0.1-20190611-2juhkwlddpdydl2yasmuia2slpev6fl5/bin/mpif90 \
	-DMPI_C_COMPILER=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-9.1.0/spectrum-mpi-10.3.0.1-20190611-2juhkwlddpdydl2yasmuia2slpev6fl5/bin/mpicc		
	# -DTPL_ARPACK_LIBRARIES="/ccs/home/liuyangz/my_software/arpack-ng-gnu/build/lib/libarpack.so;/ccs/home/liuyangz/my_software/arpack-ng-gnu/build/lib/libparpack.so"
# -DCMAKE_Fortran_FLAGS="-ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check " \
