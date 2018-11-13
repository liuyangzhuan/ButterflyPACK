#ifndef HODLR_WRAP /* allow multiple inclusions */
#define HODLR_WRAP

// #include "HODLR_config.fi"

// typedef struct { double r, i; } doublecomplex;


typedef void* F2Cptr;  // pointer passing fortran derived types to c
typedef void* C2Fptr;  // pointer passing c objects to fortran


//------------------------------------------------------------------------------
// Declartion of FORTRAN subroutines to HODLR code
extern "C" {
	
    void c_hodlr_construct(int* Npo, int* Ndim, double* Locations, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* ho_bf_for, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr), C2Fptr C_QuantApp, MPI_Fint* MPIcomm);	

    void c_hodlr_construct_matvec_init(int* Npo, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* ho_bf_for, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree);
	
    void c_hodlr_construct_matvec_compute(F2Cptr* ho_bf_for, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncHMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr), C2Fptr C_QuantApp);
	
	void c_hodlr_factor(F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree, F2Cptr* msh);	

	void c_hodlr_solve(C_DT* x, C_DT* b, int* Nloc, int* Nrhs, F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	
	
	void c_hodlr_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);	

	void c_hodlr_inv_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* ho_bf_for,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
	
	void c_hodlr_createptree(int* nmpi, int* groupmembers, MPI_Fint* MPIcomm, F2Cptr* ptree);
	
	void c_hodlr_createstats(F2Cptr* stats);		
	void c_hodlr_printstats(F2Cptr* stats,F2Cptr* ptree);		
	void c_hodlr_createoption(F2Cptr* option);	
	void c_hodlr_setoption(F2Cptr* option, char const * nam, C2Fptr val);	
	void c_hodlr_copyoption(F2Cptr* option, F2Cptr* option1);	
	
	void c_hodlr_deletestats(F2Cptr* stats);
	void c_hodlr_deleteproctree(F2Cptr* ptree);
	void c_hodlr_deletemesh(F2Cptr* msh);
	void c_hodlr_deletekernelquant(F2Cptr* ker);
	void c_hodlr_deletehobf(F2Cptr* ho_bf);
	void c_hodlr_deleteoption(F2Cptr* option);
	
	inline void c_hodlr_set_I_option(F2Cptr* option, char const * nam, int val){
		c_hodlr_setoption(option, nam, (C2Fptr) &val);
	}
	inline void c_hodlr_set_D_option(F2Cptr* option, char const * nam, double val){
		c_hodlr_setoption(option, nam, (C2Fptr) &val);
	}	


	void c_bf_construct_matvec_init(int* M, int* N,int* M_loc,int* N_loc, F2Cptr* mshr,F2Cptr* mshc,F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree);
	
	void c_bf_construct_matvec_compute(F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree,void (*C_FuncBMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr,C_DT*, C_DT*), C2Fptr C_QuantApp);	
	
	void c_bf_deletebf(F2Cptr* bf);
	
	void c_bf_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bf,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree, C_DT *a, C_DT* b);		
	
	
}
// -----------------------------------------------------------------------------



#endif
