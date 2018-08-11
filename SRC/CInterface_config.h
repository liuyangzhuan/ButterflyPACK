#ifndef __COMPLX /* allow multiple inclusions */
#define __COMPLX
typedef struct { double r, i; } doublecomplex;
#endif

#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FC_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define FC_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

typedef void* F2Cptr;  // pointer passing fortran derived types to c
typedef void* C2Fptr;  // pointer passing c objects to fortran



//------------------------------------------------------------------------------
// Declartion of FORTRAN subroutines to HODLR code
extern "C" {
	

    void FC_GLOBAL_(c_hodlr_fill,C_HODLR_FILL)(int* Npo, int* Ndim, double* Locations, int* Nmin, double* tol, int* nth, int* nmpi, int* ninc, int* preorder, int* tree, int* perms, int* Npo_loc, F2Cptr* ho_bf_for, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, doublecomplex*,C2Fptr), C2Fptr C_QuantZmn, MPI_Fint* MPIcomm);	
 
	void FC_GLOBAL_(c_hodlr_factor,C_HODLR_FACTOR)(F2Cptr* ho_bf_for,F2Cptr* ho_bf_inv,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	

	void FC_GLOBAL_(c_hodlr_solve,C_HODLR_SOLVE)(doublecomplex* x, doublecomplex* b, int* Nloc, int* Nrhs, F2Cptr* ho_bf_for,F2Cptr* ho_bf_inv,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	
	
	
	// void FC_GLOBAL_(h_matrix_apply,H_MATRIX_APPLY)(int* Npo, int* Ncol, double* Xin, double* Xout);		
}
// -----------------------------------------------------------------------------



#endif
