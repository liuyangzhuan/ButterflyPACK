This directory contains sample programs to illustrate how to call functions in ButterflyPACK from your application code. ButterflyPACK provides both double and double-complex data types. It provides both Fortran and C/C++ interfaces.

## Fortran Interface
The following pseudo codes explain how to perform construction, factorization and solve of a linear system "Z" using the Fortran interface

First, specify the data type of your application:
```
!**** DAT 1: double precision; DAT 0: double complex precision
#define DAT 1
#include "ButterflyPACK_config.fi"
```

ButterflyPACK provides two ways of constructing a hierarchical matrix.
The first option requires a user-provided function to sample any individual element of the matrix that takes the following argument list
```
subroutine Element(m,n,val,quant)
implicit none
	class(*),pointer :: quant ! quant is a user-defined derived type consisting all data and metadata needed for this user-defined function
	integer:: m,n  ! m,n specify the row and column indices of the desired matrix element
	DT::val ! val returns Z(m,n), (DT=real(kind=8) or complex(kind=8) depending on your application)

	! write your matrix element evaluation function here

end subroutine Element
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

	! write your matrix multiplication function here

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
ker%FuncHMatVec => MatVec  ! Note that at least one of ker%FuncZmn and ker%FuncHMatVec needs to set

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
call BPACK_construction_Element(bmat,option,stats,msh,ker,element_Zmn_user,ptree) ! construct a hierarchical matrix with fast matrix entry evaluation
call BPACK_construction_Matvec(bmat,matvec_user,memory,error,option,stats,ker,ptree,msh) ! construct a hierarchical matrix with fast matrix-vector multiplication
```

Factorization of the hierarchical matrix:
```
call BPACK_Factorization(bmat,option,stats,ptree,msh)
```

Solve of the hierarchical matrix:
```
call BPACK_Factorization(bmat,x,b,N_loc,nrhs,option,ptree,stats)
	! bmat is the meta-data storing the compressed and factored matrix
	! N_loc is the local number of rows/columns
	! nrhs is the number of right-hand sides
	! P is the permutation vector returned
	! b of N_loc*nrhs is the local right-hand sides
	! x of N_loc*nrhs is the local solution vector
```

## Fortran Example
You can run several examples in this directory as follows. You can modify options in the driver or from command line.

KERREG_Driver.f90:
An example for kernel ridge regression (KRR) with RBF kernel. This example constructs (with entry evaluation), factor a RBF-kernel matrix and uses it for binary classifications with UCI machine learning datasets.
```
mpirun -n nmpi ./EXAMPLE/krr
```

EMCURV_Driver.f90 and EMCURV_Module.f90:
A 2D EFIE example with several built-in geometries. This example constructs (with entry evaluation), factor the EFIE matrix and solve it with plane-wave excitations.
```
mpirun -n nmpi ./EXAMPLE/ie2d
```

EMSURF_Driver.f90 and EMSURF_Module.f90:
A 3D EFIE/CFIE example for 3D PEC surfaces. This example constructs (with entry evaluation), factor the EFIE/CFIE matrix and solve it with plane-wave excitations.
```
sh ./EM3D_DATA/preprocessor_3dmesh/run_gmsh.sh ! this preprocessor generates a few 3D example meshes using Gmsh (http://gmsh.info/)
mpirun -n nmpi ./EXAMPLE/ie3d [./EM3D_DATA/preprocessor_3dmesh/sphere_2300]
```

SMAT_Driver.f90:
An example for compressing a scattering matrix between two 3D dielectric surfaces. The scattering matrix and the coordinates of each row/column is stored in file. This example first read in the full scattering matrix, then used it as entry evaluation (if explicitflag=1) or matrix-vector multiplication (if explicitflag=0) to construct the first hierarchical matrix. Then it uses the first hierarchical matrix as matrix-vector multiplication to construct the second hierarchical matrix.
```
sh ./EM3D_DATA/FULLMAT_DATA/file_merge.sh ./EM3D_DATA/FULLMAT_DATA/Smatrix.mat ! this extract a full matrix stored as Smatrix.mat
mpirun -n nmpi ./EXAMPLE/smat
```

Frontal_Driver.f90:
An example for compressing a frontal matrix from 3D poisson equations. The frontal matrix is stored in file. This example first read in the full matrix, then used it as entry evaluation (if explicitflag=1) or matrix-vector multiplication (if explicitflag=0) to construct the first hierarchical matrix. Then it uses the first hierarchical matrix as matrix-vector multiplication to construct the second hierarchical matrix.
```
mpirun -n nmpi ./EXAMPLE/frontal
```

FULLMAT_Driver.f90:
An example for compressing a randomly generated low-rank matrix (if tst=1) or a full RBF kernel matrix stored in file (if tst=2). The example first constructs a hierarchical matrix with entry evaluation, then uses the first hierarchical matrix as matrix-vector multiplication to construct a second hierarchical matrix.
```
sh ./EM3D_DATA/FULLMAT_DATA/file_merge.sh K05N4096.csv ! this extract a full matrix stored as K05N4096.csv
mpirun -n nmpi ./EXAMPLE/full
```

## C/C++ Interface
Note that ButterflyPACK supports double and double-complex data types as two independent libraries. As such, all ButterflyPACK C++ interfaces are named with the prefix "x_". x=d for double precision and x=z for double complex precision. Take double precision for example, the caller needs to first define a class/object that can perform either matrix entry evaluation or matrix vector multiplication:
```
#include "dC_BPACK_wrapper.h"
//provide a user-defined class consisting all data and metadata needed for this matrix entry evaluation and/or matvec
class C_QuantApp {
//define your data here
};
// The entry evaluation function wrapper required by the Fortran code, val returns Z(m,n), F2Cptr is an alias of void*
inline void C_FuncZmn(int *m, int *n, double *val, F2Cptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  //call your entry evaluation function defined in C_QuantApp using Q
}
// The matvec function wrapper required by the Fortran code, see "subroutine MatVec" above for the argument list
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

Construction of the hierarchical matrix with entry evaluation:
```
d_c_bpack_construct_element(&N, &Ndim, coordinates, &nlevel, clustertree, P, &N_loc, &bmat, &option, &stats, &msh, &kerregister, &ptree, &C_FuncZmn, quant, &Fcomm);
   //N is matrix dimension
   //coordinates is a double array of size N*Ndim representing Cartesian coordinates x1(1),...,x1(Ndim),x2(1),...,x2(Ndim)....
   //if Ndim=0, coordinates is not referenced
   //clustertree is an integer array of size 2^nlevel containing leafsizes in a user-provided cluster tree
   //if nlevel=0, input requires tree(1)=N
   //P is an integer array of size N, representing permutation vector returned by the ButterflyPACK clustering
   //N_loc is the local matrix dimension
```

Construction of the hierarchical matrix with matrix-vector multiplication:
```
d_c_bpack_construct_matvec_init(&N, &nlevel, clustertree, P, &N_loc, &bmat, &option, &stats, &msh, &kerregister, &ptree);
//returning clustertree, P and N_loc can be used to define your matvec if needed
d_c_bpack_construct_matvec_compute(&bmat, &option, &stats, &msh, &kerregister, &ptree, &C_FuncHMatVec, quant);
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

In addition, ButterflyPACK can be invoked from STRUMPACK with more compact C++ interfaces than the above. See http://portal.nersc.gov/project/sparse/strumpack/ for details.

## C/C++ Example
There is one C++ example in this directory as follows. You can modify options in the driver or from command line.

InterfaceTest.cpp:
An example using a collection of simple kernels. The kernels are selected by the parameter "ker". ker=1: RBF, 2: R^4, 3: sqrt(R^2+h), 4: 1/sqrt(R^2+h), 5: (X^tY+h)^2, 6: product of two random matrices. The coordinates are generated as follows: tst=1 read from a UCI machine learning dataset, tst=2 generates a random collection of coordinates, tst=3 no coordinates are generated.
This example first constructs (with entry evaluation), factor a hierarchical matrix and then uses it as matrix-vector multiplication to construct a second hierarchical matrix.
```
mpirun -n nmpi ./EXAMPLE/ctest
```
