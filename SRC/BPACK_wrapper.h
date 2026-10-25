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
#include "bpack_gpu.h" /* the GPU evaluators of the entries (GPU_INTERFACE) */
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

void c_bpack_construct_element_compute_fortran(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bpack_construct_element_compute(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bpack_construct_init_fortran(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncDistmn)(int*, int*, double*,C2Fptr), void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);
// nproc, rank_idxs, rank_idxe: every rank's 1-based inclusive index range, so
// the mesh can carry a basis_group reproducing the H2's own distribution; a
// butterfly block built on this mesh (c_bf_construct_init) then lands on the
// same rows.
void c_bpack_set_mesh_h2(int* Npo, int* new2old, int* idxs, int* idxe, int* nproc, int* rank_idxs, int* rank_idxe, int* nmin_leaf, F2Cptr* msh_Cptr);
void c_bpack_set_mesh_distributed64(const int64_t* N_global, const int64_t* N_input_local, const int64_t* N_internal_local, const int64_t* internal_global_start, F2Cptr* msh_Cptr);
void c_bpack_wrap_h2(F2Cptr* bmat_Cptr, C2Fptr h2_ptr);
void c_bpack_get_h2(C2Fptr bmat_Cptr, F2Cptr* h2_ptr);
void c_bpack_get_gpu(C2Fptr bmat_Cptr, F2Cptr* gpu_ptr);
void c_bpack_set_gpu(C2Fptr bmat_Cptr, C2Fptr gpu_ptr);
void c_bpack_construct_init(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncDistmn)(int*, int*, double*,C2Fptr), void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);

/*
 * Distributed 64-bit ordering interface (currently implemented for format=7).
 * Global IDs are one-based int64_t values. input_locations is point-major with
 * Ndim coordinates per local point. When bounds_provided is nonzero,
 * global_bounds contains [min0,max0,...] and has length 2*Ndim.
 * c_bpack_solve/c_bpack_mult continue to accept and return vectors in the
 * caller's input rank distribution; the library redistributes them internally.
 */
typedef void (*c_bpack_func_zmn64)(const int64_t*, const int64_t*, const double*, const double*, const int*, C_DT*, C2Fptr);
typedef void (*c_bpack_func_zmn_block64)(const int64_t*, const int64_t*, const int64_t*, const int64_t*, const double*, const double*, const int*, C_DT*, const int64_t*, C2Fptr);
void c_bpack_construct_init_distributed64(const int64_t* N_global, const int64_t* N_input_local, const int* Ndim, const int64_t* input_global_ids, const double* input_locations, const int* bounds_provided, const double* global_bounds, int64_t* N_internal_local, F2Cptr* bmat, F2Cptr* option, F2Cptr* stats, F2Cptr* msh, F2Cptr* ker, F2Cptr* ptree);
void c_bpack_construct_element_compute_distributed64(F2Cptr* bmat, F2Cptr* option, F2Cptr* stats, F2Cptr* msh, F2Cptr* ker, F2Cptr* ptree, c_bpack_func_zmn64 C_FuncZmn, c_bpack_func_zmn_block64 C_FuncZmnBlock, C2Fptr C_QuantApp);
void c_bpack_get_extraction_profile(double* vals);
void c_bpack_get_distributed_layout64(F2Cptr* bmat, int64_t* N_global, int64_t* N_input_local, int64_t* N_internal_local, int64_t* internal_global_start);
/*
 * The matrix's entries on the GPU (H2 with H2_use_gpu, HODLR with
 * HODLR_use_gpu): the application's evaluator, registered after
 * c_bpack_construct_init and before the construction.  bpack_gpu.h
 * documents the three forms and their flags; doc/gpu_kernels.md is the user
 * guide.  Without one, the H2 GPU backend runs only the Schur-update pass and
 * the HODLR GPU backend evaluates the entries on the host.  Only the double
 * and double complex libraries have a GPU backend; elsewhere these are
 * ignored.
 */
void c_bpack_set_gpu_entry_launchers(F2Cptr* bmat, const bpack_gpu_entry_launchers* launchers);
void c_bpack_set_gpu_entry_source(F2Cptr* bmat, const char* source, const double* params, const int* nparams, const int* flags);
void c_bpack_set_gpu_block_evaluator(F2Cptr* bmat, bpack_gpu_block_evaluator evaluate, void* user, const int* flags);
void c_bpack_set_gpu_list_evaluator(F2Cptr* bmat, bpack_gpu_list_evaluator evaluate, void* user, const int* flags);
void c_bpack_gpu_warm_up(F2Cptr* bmat, double* seconds);
/* GPU state of formats other than H2 (SRC/BPACK_gpu_wrapper.cpp) */
void c_bpack_gpu_create(void** gpu);
void c_bpack_gpu_set_evaluator(void* gpu, const void* evaluator);
void c_bpack_get_internal_global_ids64(F2Cptr* bmat, const int64_t* internal_local_offset, const int64_t* count, int64_t* global_ids);
void c_bpack_get_input_to_internal_map64(F2Cptr* bmat, const int64_t* input_local_offset, const int64_t* count, int64_t* global_ids, int* internal_owner_ranks, int64_t* internal_global_indices);
void c_bpack_input_to_internal(F2Cptr* bmat, const int* nrhs, const C_DT* input_values, C_DT* internal_values);
void c_bpack_internal_to_input(F2Cptr* bmat, const int* nrhs, const C_DT* internal_values, C_DT* input_values);
void c_bpack_construct_init_gram(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms, int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr),  C2Fptr C_QuantApp);
void c_bpack_construct_matvec_compute(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncHMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_factor_fortran(F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*msh);
void c_bpack_factor(F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*msh);
void c_bpack_solve_fortran(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree);
void c_bpack_solve(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree);
void c_bpack_iter_usermatvec_precon(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*ker, void (*C_FuncHMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr), void (*C_FuncHPreconMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_iter_usermatvec_precon(int* Ndim, C_DT *x, C_DT *b, int*Nloc, int*Nrhs, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*ker, void (*C_FuncHMatVec_MD)(int*, char const *, int*, int*, int*, C_DT  const*,C_DT *,C2Fptr), void (*C_FuncHPreconMatVec_MD)(int*, char const *, int*, int*, int*, C_DT  const*,C_DT *,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_construct_init(int* Ns, int* Nmax, int* Ndim, double* Locations, int* perms, int* Ns_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncNearFar_MD)(int*, int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);
void c_bpack_md_construct_element_compute(int* Ndim, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn_MD)(int*, int*, int*, C_DT *,C2Fptr), void (*C_FuncZmnBlock_MD)(int*, int*, int*, int*, int64_t*, int*, int*, C_DT *, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bpack_construct_init_diaglevel_simple(int* Npo, int* Ndim, double* Locations, int* ScatterTensor, int* level, int* nlocal_active, int* local_tensor_old, int* nblock, int* block_info, int* diag_perm_size, int* diag_tensor_perm, int* mesh_range, int* route_capacity, int* local_diag_root, int* local_diag_gid, int* local_diag_pos, F2Cptr* option, F2Cptr* ptree);
void c_bpack_get_anylevel_blocks(int* level, int* nblock, int* block_info, int* mesh_range, F2Cptr* bmat, F2Cptr* option,F2Cptr* ptree,F2Cptr* msh);
void c_bpack_md_get_local_midlevel_blocks(int* Ndim, int* nc_m, int* head_array, F2Cptr* bmat, F2Cptr* option,F2Cptr* ptree,F2Cptr* msh);
void c_bpack_md_get_anylevel_blocks(int* Ndim, int* level, int* nblock, int* block_info, int* mesh_range, F2Cptr* bmat, F2Cptr* option,F2Cptr* ptree,F2Cptr* msh);
void c_bpack_md_mult(int* Ndim, char const * trans, C_DT  const * xin, C_DT * xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree,F2Cptr* msh);
void c_bpack_md_solve(int* Ndim, C_DT *x, C_DT *b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree,F2Cptr* msh);
void c_bpack_md_new2old(int* Ndim, F2Cptr* msh, int* newidx_loc, int* oldidx);
void c_bpack_md_global_new2old(int* Ndim, F2Cptr* msh, int* newidx_glo, int* oldidx);
void c_bpack_md_old2new(int* Ndim, F2Cptr* msh, int* oldidx, int* newidx);
void c_bpack_md_localindices(int* Ndim, F2Cptr* msh, int* idxs, int* Nlocal, int* Nglobal);
void c_bpack_md_vector_local2global(F2Cptr* ptree, F2Cptr* msh, int* Ndim, int* nvec, C_DT   *b, C_DT   *b_global);

void c_bpack_singleindex_to_multiindex(int* Ndim, int* dims, int* single_index_in, int* multi_index);
void c_bpack_multiindex_to_singleindex(int* Ndim, int* dims, int* single_index_in, int* multi_index);
void c_bpack_mult_fortran(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_logdet_fortran(C_DT* phase, C_RDT* logabsdet, F2Cptr* option, F2Cptr* bmat);
void c_bpack_logdet(C_DT* phase, C_RDT* logabsdet, F2Cptr* option, F2Cptr* bmat);
void c_bpack_extractelement_fortran(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps);
void c_bpack_extractelement(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps);
//------------------------------------------------------------------------------
// The same extraction from a Bmat matrix built on a subset of the caller's ranks (collective over ptree_caller, the caller's
// process tree; the source's handles bmat, option, msh, stats and ptree point to nullptr outside the subset)
void c_bpack_extractelement_bridge_fortran(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,F2Cptr* ptree_caller, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps);
void c_bpack_extractelement_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,F2Cptr* ptree_caller, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps);
//------------------------------------------------------------------------------
// Sketched blocks (format 7 with elem_extract=1): request n < Nsketch asks for Y_n = Omega_n A(rows, columns),
// sketchdim[n] x colidx[n], column major, for the rows allrows (rowidx[n] of them, stacked as for the block
// callback), the columns allcols (colidx[n]) and Omega_n (sketchdim[n] x rowidx[n]) given by its omegaidx[n]
// entries, stacked: omega_row (1-based position among the request's rows), omega_dest (1-based row of Y_n) and
// omega_val. Y_n goes to alldat_loc on the rank owning n (pgidx and pmaps as for the block callback, 1 x 1
// grids), in request order. Collective over the matrix's ranks, like the block callback.
typedef void (*c_bpack_func_zmn_sketch)(int* Nsketch, int* Nallrows, int* Nallcols, int64_t* Nallomega, int64_t* Nalldat_loc, int* allrows, int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row, int* omega_dest, double* omega_val, C_DT* alldat_loc, int* pgidx, int* Npmap, int* pmaps, C2Fptr C_QuantApp);
// An H2 (format 7) whose entries come through a collective block callback (elem_extract=1) asks for the ID
// sketches of its kernel rows through this callback (its C_QuantApp: that of the construction), set before
// c_bpack_construct_element_compute; other formats ignore it.  A nonsymmetric H2 (option%sym=0) asks for them
// only with the transposed callback set too (its ID sketches A(B, F) from the right as well), else for the
// rows' blocks.
void c_bpack_set_sketch_callback(F2Cptr* bmat, c_bpack_func_zmn_sketch C_FuncZmnSketch);
// The transposed sketched blocks, for a nonsymmetric H2: request n asks for Y_n = Omega_n A(columns, rows)^T
// (= Omega_n A^T(rows, columns)), the arguments as for the sketch callback (Omega_n over the rows).
void c_bpack_set_sketch_transpose_callback(F2Cptr* bmat, c_bpack_func_zmn_sketch C_FuncZmnSketchT);
// The sketched blocks from a matrix on a subset of the caller's ranks (as c_bpack_extractelement_bridge): a
// compression-only H2 (format 7) sketches them from its representation; other formats extract the entries and
// sketch them on their owners.
void c_bpack_extractsketch_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,F2Cptr* ptree_caller, int* Nsketch, int* Nallrows, int* Nallcols, int64_t* Nallomega, int64_t* Nalldat_loc, int* allrows, int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row, int* omega_dest, double* omega_val, C_DT* alldat_loc, int* pgidx, int* Npmap, int* pmaps);
// The same for the transposed sketches, Y_n = Omega_n A(columns, rows)^T.
void c_bpack_extractsketch_transpose_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,F2Cptr* ptree_caller, int* Nsketch, int* Nallrows, int* Nallcols, int64_t* Nallomega, int64_t* Nalldat_loc, int* allrows, int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row, int* omega_dest, double* omega_val, C_DT* alldat_loc, int* pgidx, int* Npmap, int* pmaps);
void c_bpack_inv_mult_fortran(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_inv_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
// The two above for a Fortran caller: an H2 matrix (format 7), which
// BPACK_Inv_Mult and BPACK_Mult cannot apply, with any exception reported and
// aborted here rather than thrown through the Fortran frame.
void c_bpack_h2_inv_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_h2_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bpack_set_entry_addend(F2Cptr* ker, F2Cptr* bf, F2Cptr* ptree);
void c_bpack_clear_entry_addend(F2Cptr* ker);
void c_bpack_createptree(int*nmpi, int*groupmembers, MPI_Fint*MPIcomm, F2Cptr*ptree);
void c_bpack_get_comm(F2Cptr* ptree, int* fcomm);
void c_bpack_createstats(F2Cptr*stats);
void c_bpack_new2old(F2Cptr* msh, int* newidx_loc, int* oldidx);
void c_bpack_old2new(F2Cptr* msh, int* oldidx, int* newidx);
void c_bpack_localindices(F2Cptr* msh, int* idxs, int* nlocal, int* nglobal);
void c_bpack_printstats(F2Cptr*stats, F2Cptr*ptree);
void c_bpack_printstructure(F2Cptr* bmat, int* inverse, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree);
void c_bpack_getstats(F2Cptr*stats, char const*nam, double*val_d);
void c_bpack_setstats(F2Cptr*stats, char const*nam, double*val_d);
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
void c_bf_get_block_pid(F2Cptr bf, F2Cptr ptree, int level, int index_i, int index_j, char mode, int* pid, int* myid, int* nproc);
void c_bf_get_sub_head(F2Cptr bf, F2Cptr ptree, int* is_head, MPI_Fint* comm, int* nproc);
void c_bf_get_exchange_dest(F2Cptr bf, F2Cptr ptree, int level, char mode, char collect, int nr, int nc,
                            int idx_r, int inc_r, int idx_c, int inc_c, int* pid);
void c_bf_get_mvp_vec_ranges(F2Cptr bf, F2Cptr ptree, char chara, int* idx_r, int* inc_r, int* nr,
                             int* idx_c, int* inc_c, int* nc);
void c_bf_get_all2all_dest(F2Cptr bf, F2Cptr ptree, int level, char mode_new, int nr, int nc,
                           int idx_r, int inc_r, int idx_c, int inc_c, int* pid);
void c_bf_get_all2all_range(F2Cptr bf, F2Cptr ptree, int level, char mode_new, int* idx_r, int* inc_r, int* nr,
                            int* idx_c, int* inc_c, int* nc);
void c_bf_get_block_group(F2Cptr bf, F2Cptr ptree, int level, int index_i, int index_j, char mode,
                          int* nproc, MPI_Fint* comm, int* is_head);
void c_bf_get_vec_ranges(F2Cptr bf, F2Cptr ptree, char dir, int* idx_r, int* inc_r, int* nr, int* idx_c, int* inc_c, int* nc);
void c_bf_get_layout(F2Cptr bf, int* level_butterfly, int* level_half, int* M, int* N, int* M_loc, int* N_loc, int* nlevels, int* nblocks);
void c_bf_get_blocks(F2Cptr bf, int* lvl_nr, int* lvl_nc, int* lvl_idxr, int* lvl_idxc, int* lvl_incr, int* lvl_incc, int* blk_rows, int* blk_cols, void** blk_ptr);
void c_bf_construct_init(int* M, int* N,int* M_loc,int* N_loc, int* nnsr, int* nnsc, F2Cptr* mshr,F2Cptr* mshc,F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncDistmn)(int*, int*, double*,C2Fptr), void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp);
void c_bf_construct_element_compute(F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr), void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bf_construct_matvec_compute(F2Cptr* bf, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree,void (*C_FuncBMatVec)(char const *, int*, int*, int*, C_DT const*,C_DT*,C2Fptr,C_DT*, C_DT*), C2Fptr C_QuantApp);
void c_bf_schur_deterministic(F2Cptr* bf21, F2Cptr* bf12, F2Cptr* bmat11, F2Cptr* option11, F2Cptr* stats11, F2Cptr* ptree11, F2Cptr* bfS, F2Cptr* optionS, F2Cptr* statsS, F2Cptr* ptreeS);
void c_bf_deletebf(F2Cptr*bf);
void c_bf_mult(char const * trans, C_DT const * xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bf,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree);
void c_bf_new2old_row(F2Cptr* mshr, int* newidx_loc, int* oldidx);
void c_bf_new2old_col(F2Cptr* mshc, int* newidx_loc, int* oldidx);


/* BP handles are distinct from BF/BPACK handles. Initialization mirrors BF and
 * consumes existing row/column meshes; format must be 1, 2, 3, or 5.
 * Entry/distance callbacks use signed reordered indices (+row, -column).
 * nnsr/nnsc are packed [point][knn], 1-based reordered opposite-axis indices
 * (0 denotes no neighbour); unused arrays/callbacks may be NULL.
 * Set option cpp=1. Use fresh BP stats, retain the user context, and keep
 * format unchanged between init and compute. Mult uses column-major local
 * vectors in reordered order and supports N/T/C; input and output cannot alias.
 * Delete BP with c_bp_delete; delete its output msh/ker with the BPACK helpers.
 * The input meshes, option, stats, and ptree remain caller-owned. */
void c_bp_construct_init(int* M, int* N, int* M_loc, int* N_loc, int* nnsr, int* nnsc, F2Cptr* mshr, F2Cptr* mshc, F2Cptr* bp, F2Cptr* option, F2Cptr* stats, F2Cptr* msh, F2Cptr* ker, F2Cptr* ptree, void (*C_FuncDistmn)(int*, int*, double*, C2Fptr), void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bp_construct_element_compute(F2Cptr* bp, F2Cptr* option, F2Cptr* stats, F2Cptr* msh, F2Cptr* ker, F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*, C2Fptr), void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), C2Fptr C_QuantApp);
void c_bp_mult(char const* trans, C_DT const* xin, C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, F2Cptr* bp, F2Cptr* option, F2Cptr* stats, F2Cptr* ptree);
void c_bp_delete(F2Cptr* bp);

#ifdef __cplusplus
}
#endif
//-----------------------------------------------------------------------------

#endif
