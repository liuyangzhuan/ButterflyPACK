This directory contains doxygen configuration to generate documentation of the ButterflyPACK repository.

Using CMAKE, the documentation can be generated by:  
```
mkdir build ; cd build;
cmake .. \
	-DCMAKE_Fortran_FLAGS="" \
	-DCMAKE_CXX_FLAGS="" \
	-DTPL_BLAS_LIBRARIES="${BLAS_LIB}" \
	-DTPL_LAPACK_LIBRARIES="${LAPACK_LIB}" \
	-DTPL_SCALAPACK_LIBRARIES="${SCALAPACK_LIB}" \
	-DTPL_ARPACK_LIBRARIES="${ARPACK_LIB}" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release \
	-Denable_doc=ON
make doc
```
The documentation will be generated at build/doc/html or build/doc/latex. Note that once "Denable_doc=ON" is enabled, the code cannot be compiled as the source files are inserted with doxygen markers. 


