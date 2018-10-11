#ifndef z_HODLR_WRAP /* allow multiple inclusions */
#define z_HODLR_WRAP

// #include "HODLR_config.fi"

// typedef struct { double r, i; } doublecomplex;


typedef void* F2Cptr;  // pointer passing fortran derived types to c
typedef void* C2Fptr;  // pointer passing c objects to fortran


//------------------------------------------------------------------------------
// Declartion of FORTRAN subroutines to HODLR code
extern "C" {
	
    void z_c_hodlr_construct(int* Npo, int* Ndim, double* Locations, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* ho_bf_for, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, _Complex double *,C2Fptr), C2Fptr C_QuantZmn, MPI_Fint* MPIcomm);	
 
	void z_c_hodlr_factor(F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	

	void z_c_hodlr_solve(_Complex double * x, _Complex double * b, int* Nloc, int* Nrhs, F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	
	
	void z_c_hodlr_mult(char const * trans, _Complex double * xin, _Complex double * xout, int* Nloc, int* Ncol, F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	
	
	void z_c_hodlr_createptree(int* nmpi, int* groupmembers, MPI_Fint* MPIcomm, F2Cptr* ptree);
	
	void z_c_hodlr_createstats(F2Cptr* stats);		
	void z_c_hodlr_printstats(F2Cptr* stats,F2Cptr* ptree);		
	void z_c_hodlr_createoption(F2Cptr* option);	
	void z_c_hodlr_setoption(F2Cptr* option, char const * nam, C2Fptr val);	
	
	inline void z_c_hodlr_set_I_option(F2Cptr* option, char const * nam, int val){
		z_c_hodlr_setoption(option, nam, (C2Fptr) &val);
	}
	inline void z_c_hodlr_set_D_option(F2Cptr* option, char const * nam, double val){
		z_c_hodlr_setoption(option, nam, (C2Fptr) &val);
	}		
}
// -----------------------------------------------------------------------------



#endif
