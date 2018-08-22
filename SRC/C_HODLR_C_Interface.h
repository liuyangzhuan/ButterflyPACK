#ifndef HODLR_WRAP /* allow multiple inclusions */
#define HODLR_WRAP

// #include "hodlrbf_config.h"

// typedef struct { double r, i; } doublecomplex;


typedef void* F2Cptr;  // pointer passing fortran derived types to c
typedef void* C2Fptr;  // pointer passing c objects to fortran


//------------------------------------------------------------------------------
// Declartion of FORTRAN subroutines to HODLR code
extern "C" {
	
    void c_hodlr_construct(int* Npo, int* Ndim, double* Locations, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* ho_bf_for, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, double*,C2Fptr), C2Fptr C_QuantZmn, MPI_Fint* MPIcomm);	
 
	void c_hodlr_factor(F2Cptr* ho_bf_for,F2Cptr* ho_bf_inv,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	

	void c_hodlr_solve(double* x, double* b, int* Nloc, int* Nrhs, F2Cptr* ho_bf_for,F2Cptr* ho_bf_inv,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	
	
	void c_hodlr_createptree(int* nmpi, int* groupmembers, MPI_Fint* MPIcomm, F2Cptr* ptree);
	
	void c_hodlr_createstats(F2Cptr* stats);		
	void c_hodlr_createoption(F2Cptr* option);	
	void c_hodlr_setoption(F2Cptr* option, char const * nam, C2Fptr val);	
	
	inline void set_hodlr_I_option(F2Cptr* option, char const * nam, int val){
		c_hodlr_setoption(option, nam, (C2Fptr) &val);
	}
	inline void set_hodlr_D_option(F2Cptr* option, char const * nam, double val){
		c_hodlr_setoption(option, nam, (C2Fptr) &val);
	}		
	
	// void FC_GLOBAL_(h_matrix_apply,H_MATRIX_APPLY)(int* Npo, int* Ncol, double* Xin, double* Xout);		
}
// -----------------------------------------------------------------------------



#endif
