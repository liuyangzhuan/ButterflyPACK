/*“ButterflyPACK”Copyright(c) 2018, The Regents of the University of California, through
Lawrence Berkeley National Laboratory(subject to receipt of any required approvals from the
U.S.Dept.of Energy) .All rights reserved.

If you have questions about your rights to use or distribute this software, please contact
Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

NOTICE.This Software was developed under funding from the U.S.Department of Energy and the
U.S.Government consequently retains certain rights.As such, the U.S.Government has been
granted for itself and others acting on its behalf a paid - up, nonexclusive, irrevocable
worldwide license in the Software to reproduce, distribute copies to the public, prepare
derivative works, and perform publicly and display publicly, and to permit other to do so.

Developers:Yang Liu
(Lawrence Berkeley National Lab, Computational Research Division) .
*/

#ifndef BPACK_WRAP /* allow multiple inclusions */
#define BPACK_WRAP
#include <stdint.h>
#ifdef HAVE_MPI
#include <mpi.h>
#else
#ifndef MPI_Fint
#define MPI_Fint int
#endif
#endif


// #include "ButterflyPACK_config.fi"

//typedef struct{double r, i; }doublecomplex;
typedef void*F2Cptr; //pointer passing fortran derived types to c
typedef void*C2Fptr; //pointer passing c objects to fortran


//------------------------------------------------------------------------------
//Declartion of FORTRAN subroutines to BPACK code
#ifdef __cplusplus
extern "C"{
#endif

void c_bpack_set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0);

void c_bpack_construct_element_compute(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bpack_construct_init(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncDistmn)(int*, int*, double*,C2Fptr), void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_construct_init_gram(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr),  C2Fptr C_QuantApp);
void c_bpack_construct_matvec_compute(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncHMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_factor(F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*msh);
void c_bpack_solve(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree);
void c_bpack_tfqmr_noprecon(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*ker, void (*C_FuncHMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_tfqmr_noprecon(int* Ndim, C_DT *x, C_DT *b, int*Nloc, int*Nrhs, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*ker, void (*C_FuncHMatVec_MD)(int*, char const *, int*, int*, int*, C_DT  const*,C_DT *,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_construct_init(int* Ns, int* Nmax, int* Ndim, double* Locations, int* perms, int* Ns_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncNearFar_MD)(int*, int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_construct_element_compute(int* Ndim, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn_MD)(int*, int*, int*, C_DT *,C2Fptr), void (*C_FuncZmnBlock_MD)(int*, int*, int*, int*, int64_t*, int*, int*, C_DT *, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_get_local_midlevel_blocks(int* Ndim, int* nc_m, int* head_array, F2Cptr* bmat, F2Cptr* option,F2Cptr* ptree,F2Cptr* msh);

void c_bpack_md_mult(int* Ndim, char const * trans, C_DT  const * xin, C_DT * xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree,F2Cptr* msh);
void c_bpack_md_solve(int* Ndim, C_DT *x, C_DT *b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree,F2Cptr* msh);
void c_bpack_md_new2old(int* Ndim, F2Cptr* msh, int* newidx_loc, int* oldidx);
void c_bpack_singleindex_to_multiindex(int* Ndim, int* dims, int* single_index_in, int* multi_index);
void c_bpack_multiindex_to_singleindex(int* Ndim, int* dims, int* single_index_in, int* multi_index);
void c_bpack_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_logdet(C_DT* phase, C_RDT* logabsdet, F2Cptr* option, F2Cptr* bmat);
void c_bpack_extractelement(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps);
void c_bpack_inv_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_createptree(int*nmpi, int*groupmembers, MPI_Fint*MPIcomm, F2Cptr*ptree);
void c_bpack_get_comm(F2Cptr* ptree, int* fcomm);
void c_bpack_createstats(F2Cptr*stats);
void c_bpack_new2old(F2Cptr* msh, int* newidx_loc, int* oldidx);
void c_bpack_old2new(F2Cptr* msh, int* oldidx, int* newidx);
void c_bpack_localindices(F2Cptr* msh, int* idxs, int* nlocal, int* nglobal);
void c_bpack_printstats(F2Cptr*stats, F2Cptr*ptree);
void c_bpack_printstructure(F2Cptr* bmat, int* inverse, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree);
void c_bpack_getstats(F2Cptr*stats, char const*nam, double*val_d);
void c_bpack_createoption(F2Cptr*option);
void c_bpack_setoption(F2Cptr*option, char const*nam, C2Fptr val);
void c_bpack_getoption(F2Cptr*option, char const*nam, double*val_d);
void c_bpack_readoption(F2Cptr*option, F2Cptr*ptree, int*ii);
void c_bpack_copyoption(F2Cptr*option, F2Cptr*option1);
void c_bpack_printoption(F2Cptr*option, F2Cptr*ptree);
void c_bpack_getversionnumber(int*v_major, int*v_minor, int*v_bugfix);
void c_bpack_treeindex_merged2child(int*idx_merge, int*idx_child);
void c_bpack_deletestats(F2Cptr*stats);
void c_bpack_deleteproctree(F2Cptr*ptree);
void c_bpack_deletemesh(F2Cptr*msh);
void c_bpack_md_deletemesh(int* Ndim, F2Cptr*msh);
void c_bpack_deletekernelquant(F2Cptr*ker);
void c_bpack_delete(F2Cptr*bmat);
void c_bpack_deleteoption(F2Cptr*option);
void c_bpack_vector_global2local(F2Cptr* ptree, F2Cptr* msh, int* nvec, C_DT   *b_global, C_DT   *b);
void c_bpack_vector_local2global(F2Cptr* ptree, F2Cptr* msh, int* nvec, C_DT   *b, C_DT   *b_global);

inline void c_bpack_set_I_option(F2Cptr*option, char const*nam, int val) {
c_bpack_setoption(option, nam, (C2Fptr) &val);
}
inline void c_bpack_set_D_option(F2Cptr*option, char const*nam, double val) {
c_bpack_setoption(option, nam, (C2Fptr) &val);
}

void c_bf_extractelement(F2Cptr* blocks,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps);
void c_bf_construct_init(int* M, int* N,int* M_loc,int* N_loc, int* nnsr, int* nnsc, F2Cptr* mshr,F2Cptr* mshc,F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncDistmn)(int*, int*, double*,C2Fptr), void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);
void c_bf_construct_element_compute(F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr), void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bf_construct_matvec_compute(F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree,void (*C_FuncBMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr,C_DT*, C_DT*), C2Fptr C_QuantApp);
void c_bf_deletebf(F2Cptr*bf);
void c_bf_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bf,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bf_new2old_row(F2Cptr* mshr, int* newidx_loc, int* oldidx);
void c_bf_new2old_col(F2Cptr* mshc, int* newidx_loc, int* oldidx);


#ifdef __cplusplus
}
#endif
//-----------------------------------------------------------------------------

#endif
