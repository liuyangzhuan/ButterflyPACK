module load PrgEnv-gnu
module load cpe-cuda
module load cuda
module load cmake/git-20210830

cd ..
sed -i 's/^M$//' PrecisionPreprocessing.sh
mkdir -p build
cd build
export CRAYPE_LINK_TYPE=dynamic
rm -rf CMakeCache.txt
rm -rf DartConfiguration.tcl
rm -rf CTestTestfile.cmake
rm -rf cmake_install.cmake
rm -rf CMakeFiles
cmake .. \
    -DCMAKE_INSTALL_LIBDIR=./lib \
	-DCMAKE_Fortran_FLAGS="-DMPIMODULE" \
    -DCMAKE_CXX_FLAGS="" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_COMPILER=CC \
	-DCMAKE_C_COMPILER=cc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
	-DTPL_BLAS_LIBRARIES="/opt/cray/pe/libsci/21.08.1.2/GNU/9.1/x86_64/lib/libsci_gnu_82_mpi_mp.so" \
	-DTPL_LAPACK_LIBRARIES="/opt/cray/pe/libsci/21.08.1.2/GNU/9.1/x86_64/lib/libsci_gnu_82_mpi_mp.so" \
	-DTPL_SCALAPACK_LIBRARIES="/opt/cray/pe/libsci/21.08.1.2/GNU/9.1/x86_64/lib/libsci_gnu_82_mpi_mp.so"
make ie2d    
