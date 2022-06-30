This directory contains sample programs to illustrate how to call functions in ButterflyPACK from your application code. ButterflyPACK provides single-real, double-real, single-complex and double-complex data types. It provides both Fortran and C/C++ interfaces.

## Fortran Interface to construct, factor and solve a hierarchical matrix
The following pseudo codes explain how to perform construction, factorization and solve of a linear system "Z" using the Fortran interface

First, specify the data type of your application:
```
!**** double-complex precision
#define DAT 0
#include "zButterflyPACK_config.fi"
!**** double-real precision
#define DAT 1
#include "dButterflyPACK_config.fi"
!**** single-complex precision
#define DAT 2
#include "cButterflyPACK_config.fi"
!**** single-real precision
#define DAT 3
#include "sButterflyPACK_config.fi"
```
Note that if the above header file is not included, one can still use ButterflyPACK via adding a "x_" to all derived types, functions, subroutines, and modules. x='z' for double-complex, 'd' for double-real, 'c' for single-complex, and 's' for single-complex. 

ButterflyPACK provides two ways of constructing a hierarchical matrix.
The first option requires a user-provided function to sample any individual element of the matrix that takes the following argument list
```
subroutine Element(m,n,val,quant) ! only referenced when option%elem_extract=0
implicit none
	class(*),pointer :: quant ! quant is a user-defined derived type consisting all data and metadata needed for this user-defined function
	integer:: m,n  ! m,n specify the row and column indices of the desired matrix element
	DT::val ! val returns Z(m,n), (DT=real(kind=8), complex(kind=8), real(kind=4), complex(kind=4) depending on your application)

	! write your matrix element evaluation function here

end subroutine Element
```

Sometimes it's more efficient to compute a list of matrix blocks in one shot, then the user can provide a function with the following argument list 
```
subroutine Element_ListofBlocks(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant) ! only referenced when option%elem_extract>=1
implicit none

	class(*), pointer :: quant ! quant is a user-defined derived type consisting all data and metadata needed for this user-defined function
	integer:: Ninter ! number of blocks
	integer:: allrows(:) ! 1D array containing the global row indices (in original order starting from 1 to N) stacked together
	integer:: allcols(:) ! 1D array containing the global column indices (in original order starting from 1 to N) stacked together
	DT,target:: alldat_loc(:) ! 1D array containing the local entry values defined by pmaps (in column major) stacked together
	integer:: colidx(Ninter) ! 1D array containing sizes of columns of each block
	integer:: rowidx(Ninter) ! 1D array containing sizes of rows of each block
	integer:: pgidx(Ninter)  ! 1D array containing the process group number of each block, the number starts from 1
	integer:: Npmap			 ! number of process groups
	integer:: pmaps(Npmap, 3) ! 2D array (Npmapx3) containing number of process rows, number of process columns, and starting process ID in each block

	! write your matrix element evaluation function here

end subroutine Element_ListofBlocks
```

The second option requires a user-provided function to multiply the matrix Z with an arbitrary matrix R that takes the following argument list
```
subroutine MatVec(trans,Mloc,Nloc,num_vect,R,S,quant)
implicit none
	character trans ! trans='N': S=ZR; trans='T': S^t=R^tZ; trans='C': S^c=R^cZ;
	DT:: R(:,:),S(:,:) ! represent the input matrix R and output matrix S
	integer::Mloc,Nloc ! local row and column dimensions of Z returned by the initialization subroutine BPACK_construction_Init
	integer:: num_vect ! column dimension of R
	class(*),pointer :: quant ! quant is a user-defined derived type consisting all data and metadata needed for this user-defined function

	! write your matrix multiplication function here. Note that you may need the permutation vector returned by BPACK_construction_Init to convert indices to the original order. 

end subroutine MatVec
```

Initialize ButterflyPACK metadata kernel register (ker), process tree (ptree), statistics (stats), user-options (option):
```
!*** create process tree, statistics and user-option
call CreatePtree(nmpi,groupmembers,MPI_Comm,ptree) ! groupmembers is an integer array of size nmpi, holding MPI ranks from the communicator MPI_Comm that will be used for this matrix operation
call InitStat(stats)
call SetDefaultOptions(option) ! besides the default options, other options can be set after calling this, please refer to definition of option for details

!**** register the user-defined function and derived-type in kernel register
ker%QuantApp => quant
ker%FuncZmn => Element 
ker%FuncZmnBlock => Element_ListofBlocks
ker%FuncHMatVec => MatVec  ! Note that at least one of ker%FuncZmn, Element_ListofBlocks and ker%FuncHMatVec needs to set

!**** initialization of the construction phase
call BPACK_construction_Init(N,P,N_loc,bmat,option,stats,msh,ker,ptree,Coordinates,clustertree)
	! N is matrix dimension
	! P is the permutation vector returned
	! N_loc is the local number of rows/columns
	! bmat is the meta-data storing the compressed matrix
	! Coordinates(optional) of dimension dim*N is the array of Cartesian coordinates corresponding to each row or column
	! clustertree(optional) is an array of leafsizes in a user-provided cluster tree. clustertree has length 2*nl with nl denoting level of the clustertree.
	! If clustertree is incomplete with 0 element, ButterflyPACK will adjust it to a complete tree and return a modified clustertree.
	! If the hierarchical matrix has more levels than clustertree, the code will generate more levels according to option%xyzsort, option%nogeo, and option%Nmin_leaf
```

Construction of the hierarchical matrix:
```
!**** computation of the construction phase, one the following two functions can be called depending on how ker is registered
call BPACK_construction_Element(bmat,option,stats,msh,ker,ptree) ! construct a hierarchical matrix with fast matrix entry evaluation
call BPACK_construction_Matvec(bmat,matvec_user,memory,error,option,stats,ker,ptree,msh) ! construct a hierarchical matrix with fast matrix-vector multiplication
	! matvec_user should be exactly one of z_matvec_user, d_matvec_user, c_matvec_user, s_matvec_user
	! error is the estimated error of the resulting compression
```

Factorization of the hierarchical matrix:
```
call BPACK_Factorization(bmat,option,stats,ptree,msh)
```

Solve of the hierarchical matrix:
```
call BPACK_Solution(bmat,x,b,N_loc,nrhs,option,ptree,stats)
	! bmat is the meta-data storing the compressed and factored matrix
	! N_loc is the local number of rows/columns
	! nrhs is the number of right-hand sides
	! P is the permutation vector returned
	! b of N_loc*nrhs is the local right-hand sides
	! x of N_loc*nrhs is the local solution vector
```

## Fortran Example
A number of examples are available in the build/EXAMPLE folder.
You can run several examples in this directory as follows. You can modify options in the driver or from command line using "-option --name1 val1 --name2 val2". See available command-line options from "ReadOption()". 

KERREG_Driver.f90:
(double-real) An example for kernel ridge regression (KRR) with RBF kernel. This example constructs (with entry evaluation), factor a RBF-kernel matrix and uses it for binary classifications with UCI machine learning datasets.
```
mpirun -n nmpi ./EXAMPLE/krr
```

EMCURV_Driver.f90 and EMCURV_Module.f90:
(double-complex) A 2D EFIE example with several built-in geometries. This example constructs (with entry evaluation), factor the EFIE matrix and solve it with plane-wave excitations.
```
mpirun -n nmpi ./EXAMPLE/ie2d
```

EMCURV_Eigen_Driver.f90 and EMCURV_Module.f90:
A 2D EFIE example with several built-in geometries. This example constructs (with entry evaluation), factor the EFIE matrix and compute its eigen values with ARPACK. When quant%CMmode=0, the example performs eigen analysis; when quant%CMmode=1, the example performs characteristic mode analysis. When quant%SI=0, regular mode in arpack is invoked; when quant%SI=1, shift-and-invert mode in arpack is invoked.
```
mpirun -n nmpi ./EXAMPLE/ie2deigen
```

EMSURF_Driver.f90 and EMSURF_Module.f90:
(double-complex) A 3D EFIE/CFIE example for 3D PEC surfaces. This example constructs (with entry evaluation), factor the EFIE/CFIE matrix and solve it with plane-wave excitations.
```
sh ./EM3D_DATA/preprocessor_3dmesh/run_gmsh.sh ! this preprocessor generates a few 3D example meshes using Gmsh (http://gmsh.info/)
mpirun -n nmpi ./EXAMPLE/ie3d
```

EMSURF_Driver_sp.f90 and EMSURF_Module_sp.f90:
(single-complex) A 3D EFIE/CFIE example for 3D PEC surfaces. This example constructs (with entry evaluation), factor the EFIE/CFIE matrix and solve it with plane-wave excitations.
```
sh ./EM3D_DATA/preprocessor_3dmesh/run_gmsh.sh ! this preprocessor generates a few 3D example meshes using Gmsh (http://gmsh.info/)
mpirun -n nmpi ./EXAMPLE/ie3d_sp
```

SMAT_Driver.f90:
(double-complex) An example for compressing a scattering matrix between two 3D dielectric surfaces. The scattering matrix and the coordinates of each row/column is stored in file. This example first read in the full scattering matrix, then used it as entry evaluation (if explicitflag=1) or matrix-vector multiplication (if explicitflag=0) to construct the first hierarchical matrix. Then it uses the first hierarchical matrix as matrix-vector multiplication to construct the second hierarchical matrix.
```
sh ./FULLMAT_DATA/file_merge.sh ./FULLMAT_DATA/Smatrix.mat ! this extract a full matrix stored as Smatrix.mat
mpirun -n nmpi ./EXAMPLE/smat
```

Frontal_Driver.f90:
(double-real) An example for compressing a frontal matrix from 3D poisson equations. The frontal matrix is stored in file. This example first read in the full matrix, then used it as entry evaluation (if explicitflag=1) or matrix-vector multiplication (if explicitflag=0) to construct the first hierarchical matrix. Then it uses the first hierarchical matrix as matrix-vector multiplication to construct the second hierarchical matrix.
```
mpirun -n nmpi ./EXAMPLE/frontal
```

FrontalDist_Driver.f90:
(double-complex) An example for compressing a frontal matrix from 3D elastic Helmholtz equations. The frontal matrix is stored in parallel files with a 2x2 process grid. This example first read in the 2D-block cyclic matrix, then used it as matrix-vector multiplication to construct the first hierarchical matrix. Then it uses the first hierarchical matrix as matrix-vector multiplication to construct the second hierarchical matrix.
```
sh ./FULLMAT_DATA/file_merge.sh ./FULLMAT_DATA/Frontal_elastic ! this extract 4 files storing the full matrix
mpirun -n 4 ./EXAMPLE/frontaldist
```

FULLMAT_Driver.f90:
(double-real) An example for compressing a randomly generated low-rank matrix (if tst=1) or a full RBF kernel matrix stored in file (if tst=2). The example first constructs a hierarchical matrix with entry evaluation, then uses the first hierarchical matrix as matrix-vector multiplication to construct a second hierarchical matrix.
```
sh ./FULLMAT_DATA/file_merge.sh K05N4096.csv ! this extract a full matrix stored as K05N4096.csv
mpirun -n nmpi ./EXAMPLE/full
```

## C/C++ Interface to construct, factor and solve a hierarchical matrix
Similar to the Fortran interface, the C/C++ interface is named with the prefix "x_". x='z' for double-complex, 'd' for double-real, 'c' for single-complex, and 's' for single-complex. Take double-real precision for example, the caller needs to first define a class/object that can perform either matrix entry evaluation or matrix vector multiplication:
```
#include "dC_BPACK_wrapper.h"
//provide a user-defined class consisting all data and metadata needed for this matrix entry evaluation and/or matvec
class C_QuantApp {
//define your data or metadata here
};
// The entry evaluation function wrapper required by the Fortran code, val returns Z(m,n), F2Cptr is an alias of void*. 1<=m<=N, 1<=n<=N.
inline void C_FuncZmn(int *m, int *n, double *val, F2Cptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  //call your entry evaluation function defined in C_QuantApp using Q
}

// The entry evaluation function wrapper required by the Fortran code, alldat_loc returns all local entries of Z, F2Cptr is an alias of void*.
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
	// allrows: 1D array containing the global row indices (in original order starting from 1 to N) stacked together
	// allcols: 1D array containing the global column indices (in original order starting from 1 to N) stacked together
	// alldat_loc: 1D array containing the local entry values defined by pmaps (in column major) stacked together
	// colidx: 1D array containing sizes of columns of each block
	// rowidx: 1D array containing sizes of rows of each block
	// pgidx:  1D array containing the process group number of each block, the number starts from 0
	// Npmap:  number of process groups
	// pmaps:  2D array (Npmapx3) containing number of process rows, number of process columns, and starting process ID in each block

  //call your entry evaluation function defined in C_QuantApp using Q
}

// The matvec function wrapper required by the Fortran libraries, see "subroutine MatVec" of the Fortran interface for meanings of the arguments
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  //call your matvec function defined in C_QuantApp using Q
}
```

Initialize ButterflyPACK metadata:
```
F2Cptr bmat;  //hierarchical matrix returned by Fortran code
F2Cptr option;     //option structure returned by Fortran code
F2Cptr stats;      //statistics structure returned by Fortran code
F2Cptr msh;		   //mesh structure returned by Fortran code
F2Cptr kerregister;   //kernel register quantities structure returned by Fortran code
F2Cptr ptree;      //process tree returned by Fortran code
MPI_Fint Fcomm;  // the Fortran MPI communicator
Fcomm = MPI_Comm_c2f(Comm); //Comm is the C MPI communicator provided by the user

// initialize ButterflyPACK metadata
d_c_bpack_createptree(&nmpi, groups, &Fcomm, &ptree); //groups is a int array of size nmpi holding MPI ranks out of communicator Fcomm used for ButterflyPACK
d_c_bpack_createoption(&option);
d_c_bpack_createstats(&stats);

// set ButterflyPACK options other than the default ones, the double and integer options are set with different functions:
d_c_bpack_set_D_option(&option, "name", val); //double-valued option
d_c_bpack_set_I_option(&option, "name", val); //int-valued option
```


Initialization of the construction phase
```
d_c_bpack_construct_init(&N, &Ndim, coordinates, nns_ptr, &nlevel, clustertree, P, &N_loc, &bmat, &option, &stats, &msh, &kerregister, &ptree, &C_FuncDistmn, &C_FuncNearFar, quant);
   //N is matrix dimension
   //coordinates is a double array of size N*Ndim representing Cartesian coordinates x1(1),...,x1(Ndim),x2(1),...,x2(Ndim)....
   //if Ndim=0, coordinates is not referenced
   //nns_ptr is a int array of size N*option%knn representing the knn nearest neighouring points of each of the N points. If option%knn=0, nns_ptr is not referenced.   
   //clustertree is an integer array of size 2^nlevel containing leafsizes in a user-provided cluster tree
   //if nlevel=0, input requires tree(1)=N
   //P is an integer array of size N, representing permutation vector returned by the ButterflyPACK clustering
   //N_loc is the local matrix dimension, P and N_loc can be used to define your matvec if needed
   //C_FuncDistmn: pointer to user-provided function to compute distance between any row and column of the matrix. Not referenced if option%nogeo is not 2. 
   //C_FuncNearFar: pointer to user-provided function to determine whether a block (in permuted order) is compressible or not. Not referenced if option%nogeo is not 2. 
```


Construction of the hierarchical matrix with entry evaluation:
```
d_c_bpack_construct_element_compute(&bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &C_FuncZmnBlock, quant);
   //C_FuncZmn: pointer to user-provided function to sample mn^th entry of the matrix. Used when option%elem_extract=0.
   //C_FuncZmnBlock: pointer to user-provided function to sample a list of intersections of entries of the matrix. Used when option%elem_extract>=1.
```

Construction of the hierarchical matrix with matrix-vector multiplication:
```
d_c_bpack_construct_matvec_compute(&bmat, &option, &stats, &msh, &kerregister, &ptree, &C_FuncHMatVec, quant);
   //C_FuncHMatVec: pointer to user-provided function to multiply A and A* with vectors (see above)
```
Factorization of the hierarchical matrix:
```
d_c_bpack_factor(&bmat,&option,&stats,&ptree,&msh);
```

Solve of the hierarchical matrix:
```
d_c_bpack_solve(x,b,&N_loc,&nrhs,&bmat,&option,&stats,&ptree);
	// bmat is the meta-data storing the compressed and factored matrix
	// N_loc is the local number of rows/columns
	// nrhs is the number of right-hand sides
	// b of N_loc*nrhs is the local right-hand sides (concatenation of nrhs vectors of length N_loc)
	// x of N_loc*nrhs is the local solution vector (concatenation of nrhs vectors of length N_loc)
```
# Interface to construct, and multiply a low-rank/butterfly block
ButterflyPACK also provides an interface for constructing a single MxN low-rank/butterfly block. Take double-real precision for example, the caller needs to first define a class/object that can perform either matrix entry evaluation or matrix vector multiplication:
```
#include "dC_BPACK_wrapper.h"
//provide a user-defined class consisting all data and metadata needed for this matrix entry evaluation and/or matvec
class C_QuantApp {
//define your data or metadata here
};
// The entry evaluation function wrapper required by the Fortran code, val returns Z(m,n), F2Cptr is an alias of void*. If n<0 (column index) and m>0 (row index), then 1<=m<=M, -N<=n<=-1; if m<0 (column index) and n>0 (row index), then 1<=n<=M, -N<=m<=-1. 
inline void C_FuncBZmn(int *m, int *n, double *val, F2Cptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  //call your entry evaluation function defined in C_QuantApp using Q
}

// The entry evaluation function wrapper required by the Fortran code, alldat_loc returns all local entries of Z, F2Cptr is an alias of void*.
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
	// allrows: 1D array containing the global row indices (in original order starting from 1 to N) stacked together
	// allcols: 1D array containing the global column indices (in original order starting from 1 to N) stacked together
	// alldat_loc: 1D array containing the local entry values defined by pmaps (in column major) stacked together
	// colidx: 1D array containing sizes of columns of each block
	// rowidx: 1D array containing sizes of rows of each block
	// pgidx:  1D array containing the process group number of each block, the number starts from 0
	// Npmap:  number of process groups
	// pmaps:  2D array (Npmapx3) containing number of process rows, number of process columns, and starting process ID in each block

  //call your entry evaluation function defined in C_QuantApp using Q
}

// The matvec function wrapper required by the Fortran libraries.
inline void C_FuncBMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant, double *a, double *b){
  C_QuantApp* Q = (C_QuantApp*) quant;
  //call your matvec function defined in C_QuantApp using Q
}
```

Initialize ButterflyPACK metadata:
```
F2Cptr bf_b;  //hierarchical matrix returned by Fortran code
F2Cptr option_b;     //option structure returned by Fortran code
F2Cptr stats_b;      //statistics structure returned by Fortran code
F2Cptr msh_b;		   //mesh structure returned by Fortran code
F2Cptr kerregister_b;   //kernel register quantities structure returned by Fortran code
F2Cptr ptree_b;      //process tree returned by Fortran code
MPI_Fint Fcomm;  // the Fortran MPI communicator
Fcomm = MPI_Comm_c2f(Comm); //Comm is the C MPI communicator provided by the user

// initialize ButterflyPACK metadata
d_c_bpack_createptree(&nmpi, groups, &Fcomm, &ptree_b); //groups is a int array of size nmpi holding MPI ranks out of communicator Fcomm used for ButterflyPACK
d_c_bpack_createoption(&option_b);
d_c_bpack_createstats(&stats_b);

// set ButterflyPACK options other than the default ones, the double and integer options are set with different functions:
d_c_bpack_set_D_option(&option_b, "name", val); //double-valued option
d_c_bpack_set_I_option(&option_b, "name", val); //int-valued option
```

Initialization of the construction phase. In this interface, ButterflyPACK treats the single block as the 1-2 off-diagonal block of a larger hierarchical matrix. 
```
// one needs to first a mesh metadata for the row and column dimension, by calling d_c_bpack_construct_init. The *_dummy argument is not referenced. 
d_c_bpack_construct_init(&M, &Ndim, dat_ptr_m, nns_ptr_m,&nlevel_m, tree_m, perms_m, &myseg_m, &bmat_dummy, &option_b, &stats_dummy, &msh0_m, &kerquant_dummy, &ptree_b, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);
d_c_bpack_construct_init(&N, &Ndim, dat_ptr_n, nns_ptr_n,&nlevel_n, tree_n, perms_n, &myseg_n, &bmat_dummy, &option_b, &stats_dummy, &msh0_n, &kerquant_dummy, &ptree_b, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);
// the construction initialization can be done by:
d_c_bf_construct_init(&M, &N, &myseg_m, &myseg_n, nns_ptr_m, nns_ptr_n, &msh0_m, &msh0_n, &bf_b, &option_b, &stats_b, &msh_b, &kerregister_b, &ptree_b,&C_FuncDistmn, &C_FuncNearFar, quant);
```

Construction of the block with entry evaluation:
```
d_c_bf_construct_element_compute(&bf_b, &option_b, &stats_b, &msh_b, &kerregister_b, &ptree_b, &C_FuncBZmn, &C_FuncBZmnBlock, quant); 
//C_FuncBZmn: pointer to user-provided function to sample mn^th entry of the matrix. Used when option%elem_extract=0.
//C_FuncBZmnBlock: pointer to user-provided function to sample a list of intersections of entries of the matrix. Used when option%elem_extract>=1.
```

Construction of the block with matrix-vector multiplication:
```
d_c_bf_construct_matvec_compute(&bf_b, &option_b, &stats_b, &msh_b, &kerregister_b, &ptree_b, &C_FuncBMatVec, quant);
   //C_FuncBMatVec: pointer to user-provided function to multiply A and A* with vectors (see above)
```

Application of the constructed block:
```
c_bf_mult(trans, xin, xout, Ninloc, Noutloc, Nvec, bf_b, option_b, stats_b, ptree_b);
   // trans: 'N', 'C', or 'T'
   // xin: input vector of Ninloc*Nvec (concatenation of Nvec vectors of length Ninloc)
   // xout: output vector of Noutloc*Nvec (concatenation of Nvec vectors of length Noutloc)
   // Ninloc: length of local input vectors
   // Noutloc: length of local output vectors
   // Nvec: number of vectors
```
In addition, ButterflyPACK can be invoked from STRUMPACK with more compact C++ interfaces than the above. See http://portal.nersc.gov/project/sparse/strumpack/ for details.

## C/C++ Example
A number of examples are available in the build/EXAMPLE folder.
You can modify options in the driver or from command line.

InterfaceTest.cpp:
An example using a collection of simple kernels. The kernels are selected by the parameter "ker". ker=1: RBF, 2: R^4, 3: sqrt(R^2+h), 4: 1/sqrt(R^2+h), 5: (X^tY+h)^2, 6: product of two random matrices. The coordinates are generated as follows: tst=1 read from a UCI machine learning dataset, tst=2 generates a random collection of coordinates, tst=3 or 4: no coordinates are generated.
This example first constructs (with entry evaluation), factor a hierarchical matrix and then uses it as matrix-vector multiplication to construct a second hierarchical matrix.
```
mpirun -n nmpi ./EXAMPLE/ctest
```

Taylor2D.cpp:
An example that solves a 2D EFIE assuming inhomogeneous media. A approximate Green's function is used to evaluate any entry of the matrix. 
```
mpirun -n nmpi ./EXAMPLE/go2d
```

FIO_Driver.cpp:
An example that constructs two 1D Fourier integral operators (FIOs) as two butterfly-compressed blocks using entry-evaluation-based APIs, then construct their product as another butterfly-compressed block using matvec-based APIs. 
```
mpirun -n nmpi ./EXAMPLE/cfio
```

InverseFIO_Driver.cpp:
An example that constructs a 1D FIO as a butterlfy block A, computes a rhs by b=A*x_ref, then generates an approximate solution using (A^*A)^-1A^*b. Then matrix A^*A is compressed as an HODLR matrix using matvec-based APIs. 
```
mpirun -n nmpi ./EXAMPLE/cifio
```
*/
