/*
 * STRUMPACK -- STRUctured Matrices PACKage, Copyright (c) 2014, The
 * Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals
 * from the U.S. Dept. of Energy).  All rights reserved.
 *
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov.
 *
 * NOTICE. This software is owned by the U.S. Department of Energy. As
 * such, the U.S. Government has been granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable,
 * worldwide license in the Software to reproduce, prepare derivative
 * works, and perform publicly and display publicly.  Beginning five
 * (5) years after the date permission to assert copyright is obtained
 * from the U.S. Department of Energy, and subject to any subsequent
 * five (5) year renewals, the U.S. Government is granted for itself
 * and others acting on its behalf a paid-up, nonexclusive,
 * irrevocable, worldwide license in the Software to reproduce,
 * prepare derivative works, distribute copies to the public, perform
 * publicly and display publicly, and to permit others to do so.
 *
 * Developers: Pieter Ghysels, Francois-Henry Rouet, Xiaoye S. Li.
 *             (Lawrence Berkeley National Lab, Computational Research
 *             Division).
 *
 */
/*! \file CPPWrapper.hpp
 * \brief Templated CPP Interface to ButterflyPACK, modified from HODLRWrapper.hpp of STRUMPACK
 */
#ifndef BUTTERFLYPACK_WRAPPER_HPP
#define BUTTERFLYPACK_WRAPPER_HPP

#include <cassert>
#include <complex>
#include <cstdint>
#include <iostream>
#include <vector>
#include <memory>
#include <cstring>
#include <getopt.h>
#include <unistd.h>

typedef void* F2Cptr;
typedef void* C2Fptr;

namespace butterflypack {

  template<class T> struct RealType { typedef T value_type; };
  template<class T> struct RealType<std::complex<T>> { typedef T value_type; };

    void bpack_createptree
    (int& P, int* groups, MPI_Fint comm, F2Cptr& ptree);

    void bpack_createoptions(F2Cptr& options);

    void bpack_copyoptions(F2Cptr& in, F2Cptr& out);

    void bpack_printoptions(F2Cptr& options, F2Cptr& ptree);

    void bpack_printstats(F2Cptr& stats, F2Cptr& ptree);

    void bpack_createstats(F2Cptr& stats);

    void bpack_set_option
    (F2Cptr& options, const std::string& opt, double v);

    void bpack_set_option
    (F2Cptr& options, const std::string& opt, int v);

    /**
     * Possible values:
     *
     *  Time_Fill, Time_Factor, Time_Solve, Time_Sblock, Time_Inv,
     *  Time_SMW, Time_RedistB, Time_RedistV, Time_C_Mult,
     *  Time_Direct_LU, Time_Add_Multiply, Time_Multiply, Time_XLUM,
     *  Time_Split, Time_Comm, Time_Idle
     *
     *  Flop_Fill, Flop_Factor, Flop_Solve, Flop_C_Mult
     *
     *  Mem_Factor, Mem_Fill, Mem_Sblock, Mem_SMW, Mem_Direct_inv,
     *  Mem_Direct_for, Mem_int_vec, Mem_Comp_for
     *
     *  Rank_max
     */
    double bpack_get_stat
    (F2Cptr& stats, const std::string& name);

    template<typename scalar_t>
    void bpack_construct_ho_init
    (int N, int d, double* data, int* nns, int& lvls, int* tree, int* perm,
     int& lrow, F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata);

    template<typename scalar_t>
    void bpack_construct_ho_init_Gram
    (int N, int d, double* data, int* nns, int& lvls, int* tree, int* perm,
     int& lrow, F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, scalar_t*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, scalar_t* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata);

    template<typename scalar_t> void bpack_construct_ho_element_compute
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, scalar_t*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, scalar_t* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr K);

    template<typename scalar_t> void bpack_construct_ho_matvec_compute
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*matvec)
     (char const*, int*, int*, int*, const scalar_t*, scalar_t*, C2Fptr),
     C2Fptr fdata);

    template<typename scalar_t>
    void bpack_construct_bf_init
    (int M, int N, int& lrows, int& lcols, int* nsr, int* nnsc,
     F2Cptr rmsh, F2Cptr cmsh, F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata);

    template<typename scalar_t> void bpack_construct_bf_matvec_compute
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (const char*, int*, int*, int*, const scalar_t*,
      scalar_t*, C2Fptr, scalar_t*, scalar_t*), C2Fptr fdata);

    template<typename scalar_t> void bpack_construct_bf_element_compute
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*element)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, scalar_t* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata);

    template<typename scalar_t> void bpack_extract_elements_ho
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols,
     std::int64_t Nalldat_loc, int* allrows, int* allcols, scalar_t* alldat_loc,
     int* rowidx, int* colidx, int* pgidx, int Npmap, int* pmaps);

    template<typename scalar_t> void bpack_extract_elements_bf
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols,
     std::int64_t Nalldat_loc, int* allrows, int* allcols, scalar_t* alldat_loc,
     int* rowidx, int* colidx, int* pgidx, int Npmap, int* pmaps);

    void bpack_deletestats(F2Cptr&);

    void bpack_deleteproctree(F2Cptr&);

    void bpack_deletemesh(F2Cptr&);

    void bpack_deletekernelquant(F2Cptr&);

    void bpack_delete_ho(F2Cptr&);

    void bpack_delete_bf(F2Cptr&);

    void bpack_deleteoptions(F2Cptr&);

    void bpack_getversionnumber(int& v_major, int& v_minor, int& v_bugfix);

    template<typename scalar_t> void bpack_mult_ho
    (char op, const scalar_t* X, scalar_t* Y, int Xlrows, int Ylrows, int cols,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree);

    template<typename scalar_t> void bpack_mult_bf
    (char op, const scalar_t* X, scalar_t* Y, int Xlrows, int Ylrows, int cols,
     F2Cptr lr_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree);

    template<typename scalar_t> void bpack_factor_ho
    (F2Cptr& ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree, F2Cptr msh);

    template<typename scalar_t> void bpack_solve_ho
    (scalar_t* X, const scalar_t* B, int lrows, int rhs,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree);

    template<typename scalar_t> void bpack_inv_mult_ho
    (char op, const scalar_t* B, scalar_t* X, int Xlrows, int Blrows, int rhs,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree);

    int bpack_treeindex_merged2child(int idx_merge);
    int bpack_new2old(F2Cptr msh, int i_new_loc);
    void bpack_set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0);
} // end namespace butterflypack

#endif // BUTTERFLYPACK_WRAPPER_HPP
