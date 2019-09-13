 export TRAVIS_BUILD_DIR=$PWD
 export BLUE="\033[34;1m"
 export CXX="g++-4.8"
 export CC="gcc-4.8"

 export BLAS_LIB=/usr/lib/x86_64-linux-gnu/libblas.so
 export LAPACK_LIB=/usr/lib/x86_64-linux-gnu/liblapack.so
 export SCALAPACK_LIB=$TRAVIS_BUILD_DIR/installDir/scalapack-2.0.2/libscalapack.a

 export ARPACK_LIB="$TRAVIS_BUILD_DIR/installDir/arpack-ng/build/lib/libarpack.so;$TRAVIS_BUILD_DIR/installDir/arpack-ng/build/lib/libparpack.so"
 export BLUE="\033[34;1m"
 rm -rf build
 mkdir -p build
 cd build
 rm -rf CMakeCache.txt
 rm -rf DartConfiguration.tcl
 rm -rf CTestTestfile.cmake
 rm -rf cmake_install.cmake
 rm -rf CMakeFiles
    cmake .. \
    -DCMAKE_CXX_FLAGS="-fopenmp" \
    -Denable_blaslib=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_C_COMPILER=/usr/bin/mpicc \
    -DCMAKE_CXX_COMPILER=/usr/bin/mpic++ \
    -DCMAKE_Fortran_COMPILER=/usr/bin/mpif90 \
    -DCMAKE_INSTALL_PREFIX=. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DCMAKE_Fortran_FLAGS="-fopenmp" \
    -DTPL_BLAS_LIBRARIES="$BLAS_LIB" \
    -DTPL_LAPACK_LIBRARIES="$LAPACK_LIB" \
    -DTPL_SCALAPACK_LIBRARIES="$SCALAPACK_LIB"
 #   -DTPL_ARPACK_LIBRARIES="$ARPACK_LIB"
 make
 make install


