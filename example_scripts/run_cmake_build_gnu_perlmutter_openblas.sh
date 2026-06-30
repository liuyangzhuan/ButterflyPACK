module load PrgEnv-gnu
module load cmake
module load python/3.11
module load cray-fftw 
module load cray-libsci

cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
export ZFP_INSTALL_DIR=$HOME/zfp-install
#export PARSEC_INSTALL_DIR=/global/cfs/cdirs/m2957/liuyangz/my_software/parsec_pr759/install
#export CMAKE_PREFIX_PATH="$PARSEC_INSTALL_DIR:${CMAKE_PREFIX_PATH}"
#export PATH="$PARSEC_INSTALL_DIR/bin:${PATH}"
#export LD_LIBRARY_PATH="$PARSEC_INSTALL_DIR/lib64:${LD_LIBRARY_PATH}"

rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
rm -rf SRC_DOUBLE SRC_DOUBLECOMPLEX SRC_SINGLE SRC_COMPLEX
cmake .. \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE" \
    -DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-Denable_python=OFF \
	-Denable_parsec=OFF \
	-Denable_toplevel_openmp=OFF \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DOpenMP_C_FLAGS="-fopenmp" \
	-DOpenMP_C_LIB_NAMES="gomp" \
	-DOpenMP_Fortran_FLAGS="-fopenmp" \
	-DOpenMP_Fortran_LIB_NAMES="gomp" \
	-DOpenMP_omp_LIBRARY=$(gcc --print-file-name=libgomp.so) \
	-DTPL_ZFP_LIBRARIES="$ZFP_INSTALL_DIR/lib64/libzFORp.so;$ZFP_INSTALL_DIR/lib64/libzfp.so" \
	-DTPL_ZFP_INCLUDE="$ZFP_INSTALL_DIR/include" \
	-DTPL_BLAS_LIBRARIES="$CRAY_LIBSCI_PREFIX_DIR/lib/libsci_gnu_mp.so" \
	-DTPL_LAPACK_LIBRARIES="$CRAY_LIBSCI_PREFIX_DIR/lib/libsci_gnu_mp.so" \
	-DTPL_SCALAPACK_LIBRARIES="$CRAY_LIBSCI_PREFIX_DIR/lib/libsci_gnu_mp.so"

make ie2d -j16
make ie3d -j16
make ie3d_sp -j16
make ie3d_mp -j16
make ie3dport -j16
make frankben -j16   
make frankben_t -j16
make cifio2dsb -j16  
make cvie2d -j16  
make cvie2d_t -j16


#       -DPaRSEC_DIR="$PARSEC_INSTALL_DIR/share/cmake/parsec" \
