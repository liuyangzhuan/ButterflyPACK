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
/*! \file CPPWrapper.cpp
 * \brief Templated CPP Interface to ButterflyPACK, modified from HODLRWrapper.cpp of STRUMPACK
 */
#include <cassert>
#include <complex>
#include <iomanip>

#define OMPI_SKIP_MPICXX 1
#ifdef HAVE_MPI
#include <mpi.h>
#else
#ifndef MPI_Fint
#define MPI_Fint int
#endif
#endif

#include "CPPWrapper.hpp"
#include "sC_BPACK_wrapper.h"
#include "dC_BPACK_wrapper.h"
#include "cC_BPACK_wrapper.h"
#include "zC_BPACK_wrapper.h"

namespace butterflypack {

    void bpack_createptree
    (int& P, int* groups, MPI_Fint comm, F2Cptr& ptree) {
      d_c_bpack_createptree(&P, groups, &comm, &ptree);
    }

    void bpack_createoptions(F2Cptr& options) {
      d_c_bpack_createoption(&options);
    }


    void bpack_copyoptions(F2Cptr& in, F2Cptr& out) {
      d_c_bpack_copyoption(&in, &out);
    }

    void bpack_printoptions(F2Cptr& options, F2Cptr& ptree) {
      d_c_bpack_printoption(&options, &ptree);
    }

    void bpack_printstats(F2Cptr& stats, F2Cptr& ptree) {
      d_c_bpack_printstats(&stats, &ptree);
    }

    void bpack_createstats(F2Cptr& stats) {
      d_c_bpack_createstats(&stats);
    }

    void bpack_set_option(F2Cptr& options, const std::string& opt, double v) {
      d_c_bpack_set_D_option(&options, opt.c_str(), v);
    }

    void bpack_set_option(F2Cptr& options, const std::string& opt, int v) {
      d_c_bpack_set_I_option(&options, opt.c_str(), v);
    }


    double bpack_get_stat
    (F2Cptr stats, const std::string& name) {
      double val;
      d_c_bpack_getstats(&stats, name.c_str(), &val);
      return val;
    }


    template<> void bpack_construct_ho_init<float>
    (int N, int d, double* data, int* nns, int lvls, int* tree, int* perm,
     int& lrow, F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
        s_c_bpack_construct_init
          (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
           &stats, &msh, &kerquant, &ptree, C_FuncDistmn, C_FuncNearFar,
           fdata);
    }
    template<> void bpack_construct_ho_init<double>
    (int N, int d, double* data, int* nns, int lvls, int* tree, int* perm,
     int& lrow, F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      d_c_bpack_construct_init
        (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncDistmn, C_FuncNearFar,
         fdata);
    }
    template<> void bpack_construct_ho_init<_Complex float>
    (int N, int d, double* data, int* nns, int lvls, int* tree,
     int* perm, int& lrow, F2Cptr& ho_bf, F2Cptr& options,
     F2Cptr& stats, F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      c_c_bpack_construct_init
        (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
          &stats, &msh, &kerquant, &ptree, C_FuncDistmn,
          C_FuncNearFar, fdata);
    }
    template<> void bpack_construct_ho_init<_Complex double>
    (int N, int d, double* data, int* nns, int lvls, int* tree,
     int* perm, int& lrow, F2Cptr& ho_bf, F2Cptr& options,
     F2Cptr& stats, F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      z_c_bpack_construct_init
        (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncDistmn,
         C_FuncNearFar, fdata);
    }

    template<> void bpack_construct_ho_init_Gram<float>
    (int N, int d, double* data, int* nns, int lvls, int* tree, int* perm,
     int& lrow, F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, float*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, float* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
        s_c_bpack_construct_init_gram
          (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
           &stats, &msh, &kerquant, &ptree, C_FuncZmn, C_FuncZmnBlock, fdata);
    }
    template<> void bpack_construct_ho_init_Gram<double>
    (int N, int d, double* data, int* nns, int lvls, int* tree, int* perm,
     int& lrow, F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, double* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      d_c_bpack_construct_init_gram
        (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncZmn, C_FuncZmnBlock, fdata);
    }
    template<> void bpack_construct_ho_init_Gram<_Complex float>
    (int N, int d, double* data, int* nns, int lvls, int* tree,
     int* perm, int& lrow, F2Cptr& ho_bf, F2Cptr& options,
     F2Cptr& stats, F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, _Complex float*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, _Complex float* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
        c_c_bpack_construct_init_gram
          (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
           &stats, &msh, &kerquant, &ptree,
           reinterpret_cast<
           void(*)(int*, int*, _Complex float*, C2Fptr)>(C_FuncZmn),
           reinterpret_cast<
           void(*)(int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
                   int* allrows, int* allcols, _Complex float* alldat_loc,
                   int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
                   C2Fptr elems)>(C_FuncZmnBlock), fdata);
    }
    template<> void bpack_construct_ho_init_Gram<_Complex double>
    (int N, int d, double* data, int* nns, int lvls, int* tree,
     int* perm, int& lrow, F2Cptr& ho_bf, F2Cptr& options,
     F2Cptr& stats, F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, _Complex double*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, _Complex double* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      z_c_bpack_construct_init_gram
        (&N, &d, data, nns, &lvls, tree, perm, &lrow, &ho_bf, &options,
         &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<
         void(*)(int*, int*, _Complex double*, C2Fptr)>(C_FuncZmn),
         reinterpret_cast<
         void(*)(int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
                 int* allrows, int* allcols, _Complex double* alldat_loc,
                 int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
                 C2Fptr elems)>(C_FuncZmnBlock), fdata);
    }


    template<> void bpack_construct_ho_element_compute<float>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, float*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, float* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      s_c_bpack_construct_element_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree,
         C_FuncZmn, C_FuncZmnBlock, fdata);
    }
    template<> void bpack_construct_ho_element_compute<double>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, double* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      d_c_bpack_construct_element_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree,
         C_FuncZmn, C_FuncZmnBlock, fdata);
    }
    template<> void bpack_construct_ho_element_compute<_Complex float>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, _Complex float*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, _Complex float* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      c_c_bpack_construct_element_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<
         void(*)(int*, int*, _Complex float*, C2Fptr)>(C_FuncZmn),
         reinterpret_cast<
         void(*)(int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
                 int* allrows, int* allcols, _Complex float* alldat_loc,
                 int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
                 C2Fptr elems)>(C_FuncZmnBlock),
         fdata);
    }
    template<> void bpack_construct_ho_element_compute<_Complex double>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncZmn)(int*, int*, _Complex double*, C2Fptr),
     void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, _Complex double* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      z_c_bpack_construct_element_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<
         void(*)(int*, int*, _Complex double*, C2Fptr)>(C_FuncZmn),
         reinterpret_cast<
         void(*)(int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
                 int* allrows, int* allcols, _Complex double* alldat_loc,
                 int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
                 C2Fptr elems)>(C_FuncZmnBlock),
         fdata);
    }

    template<> void bpack_construct_ho_matvec_compute<float>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (char const*, int*, int*, int*, const float*, float*, C2Fptr),
     C2Fptr fdata) {
      s_c_bpack_construct_matvec_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree, matvec, fdata);
    }
    template<> void bpack_construct_ho_matvec_compute<double>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (char const*, int*, int*, int*, const double*, double*, C2Fptr),
     C2Fptr fdata) {
      d_c_bpack_construct_matvec_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree, matvec, fdata);
    }
    template<> void bpack_construct_ho_matvec_compute<_Complex float>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (char const*, int*, int*, int*, const _Complex float*,
      _Complex float*, C2Fptr), C2Fptr fdata) {
      c_c_bpack_construct_matvec_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<
         void(*)(char const*, int*, int*, int*, const _Complex float*,
                 _Complex float*, C2Fptr)>(matvec), fdata);
    }
    template<> void bpack_construct_ho_matvec_compute<_Complex double>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (char const*, int*, int*, int*, const _Complex double*,
      _Complex double*, C2Fptr), C2Fptr fdata) {
      z_c_bpack_construct_matvec_compute
        (&ho_bf, &options, &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<
         void(*)(char const*, int*, int*, int*, const _Complex double*,
                 _Complex double*, C2Fptr)>(matvec), fdata);
    }

    template<> void bpack_construct_bf_init<float>
    (int M, int N, int& lrows, int& lcols, int* nnsr, int* nnsc,
     F2Cptr rmsh, F2Cptr cmsh, F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      s_c_bf_construct_init
        (&M, &N, &lrows, &lcols, nnsr, nnsc, &rmsh, &cmsh, &lr_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncDistmn, C_FuncNearFar,
         fdata);
    }
    template<> void bpack_construct_bf_init<double>
    (int M, int N, int& lrows, int& lcols, int* nnsr, int* nnsc,
     F2Cptr rmsh, F2Cptr cmsh, F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      d_c_bf_construct_init
        (&M, &N, &lrows, &lcols, nnsr, nnsc, &rmsh, &cmsh, &lr_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncDistmn, C_FuncNearFar,
         fdata);
    }
    template<> void bpack_construct_bf_init<_Complex float>
    (int M, int N, int& lrows, int& lcols, int* nnsr, int* nnsc,
     F2Cptr rmsh, F2Cptr cmsh, F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      c_c_bf_construct_init
        (&M, &N, &lrows, &lcols, nnsr, nnsc, &rmsh, &cmsh, &lr_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncDistmn, C_FuncNearFar,
         fdata);
    }
    template<> void bpack_construct_bf_init<_Complex double>
    (int M, int N, int& lrows, int& lcols, int* nnsr, int* nnsc,
     F2Cptr rmsh, F2Cptr cmsh, F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats,
     F2Cptr& msh, F2Cptr& kerquant, F2Cptr& ptree,
     void (*C_FuncDistmn)(int*, int*, double*, C2Fptr),
     void (*C_FuncNearFar)(int*, int*, int*, C2Fptr), C2Fptr fdata) {
      z_c_bf_construct_init
        (&M, &N, &lrows, &lcols, nnsr, nnsc, &rmsh, &cmsh, &lr_bf, &options,
         &stats, &msh, &kerquant, &ptree, C_FuncDistmn, C_FuncNearFar,
         fdata);
    }

    template<> void bpack_construct_bf_matvec_compute<float>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (const char*, int*, int*, int*, const float*,
      float*, C2Fptr, float*, float*), C2Fptr fdata) {
      s_c_bf_construct_matvec_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         matvec, fdata);
    }
    template<> void bpack_construct_bf_matvec_compute<double>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (const char*, int*, int*, int*, const double*,
      double*, C2Fptr, double*, double*), C2Fptr fdata) {
      d_c_bf_construct_matvec_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         matvec, fdata);
    }
    template<> void bpack_construct_bf_matvec_compute<_Complex float>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (char const*, int*, int*, int*, const _Complex float*,
      _Complex float*, C2Fptr, _Complex float*,
      _Complex float*), C2Fptr fdata) {
      c_c_bf_construct_matvec_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<void(*)
         (char const*, int*, int*, int*, const _Complex float*,
          _Complex float*, C2Fptr, _Complex float*,
          _Complex float*)>(matvec), fdata);
    }
    template<> void bpack_construct_bf_matvec_compute<_Complex double>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*matvec)
     (char const*, int*, int*, int*, const _Complex double*,
      _Complex double*, C2Fptr, _Complex double*,
      _Complex double*), C2Fptr fdata) {
      z_c_bf_construct_matvec_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         reinterpret_cast<void(*)
         (char const*, int*, int*, int*, const _Complex double*,
          _Complex double*, C2Fptr, _Complex double*,
          _Complex double*)>(matvec), fdata);
    }

    template<> void bpack_construct_bf_element_compute<float>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, float* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      s_c_bf_construct_element_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         nullptr, C_FuncZmnBlock, fdata);
    }
    template<> void bpack_construct_bf_element_compute<double>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, double* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      d_c_bf_construct_element_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         nullptr, C_FuncZmnBlock, fdata);
    }
    template<> void bpack_construct_bf_element_compute<_Complex float>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, _Complex float* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      c_c_bf_construct_element_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         nullptr, reinterpret_cast<void(*)
         (int*, int*, int*, std::int64_t*, int*, int*, _Complex float*,
          int*, int*, int*, int*, int*, C2Fptr)>(C_FuncZmnBlock), fdata);
    }
    template<> void bpack_construct_bf_element_compute<_Complex double>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& stats, F2Cptr& msh,
     F2Cptr& kerquant, F2Cptr& ptree, void (*C_FuncZmnBlock)
     (int* Ninter, int* Nallrows, int* Nallcols, std::int64_t* Nalldat_loc,
      int* allrows, int* allcols, _Complex double* alldat_loc,
      int* rowids, int* colids, int* pgids, int* Npmap, int* pmaps,
      C2Fptr elems), C2Fptr fdata) {
      z_c_bf_construct_element_compute
        (&lr_bf, &options, &stats, &msh, &kerquant, &ptree,
         nullptr, reinterpret_cast<void(*)
         (int*, int*, int*, std::int64_t*, int*, int*, _Complex double*,
          int*, int*, int*, int*, int*, C2Fptr)>(C_FuncZmnBlock), fdata);
    }

    template<> void bpack_extract_elements_ho<float>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, float* alldat_loc, int* rowidx, int* colidx,
     int* pgidx, int Npmap, int* pmaps) {
      s_c_bpack_extractelement
        (&ho_bf, &options, &msh, &stats, &ptree, &Ninter, &Nallrows,
         &Nallcols, &Nalldat_loc, allrows, allcols, alldat_loc,
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<float>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }
    template<> void bpack_extract_elements_ho<double>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, double* alldat_loc, int* rowidx, int* colidx,
     int* pgidx, int Npmap, int* pmaps) {
      d_c_bpack_extractelement
        (&ho_bf, &options, &msh, &stats, &ptree, &Ninter, &Nallrows,
         &Nallcols, &Nalldat_loc, allrows, allcols, alldat_loc,
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<double>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }
    template<> void bpack_extract_elements_ho<_Complex float>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, _Complex float* alldat_loc,
     int* rowidx, int* colidx, int* pgidx, int Npmap, int* pmaps) {
      c_c_bpack_extractelement
        (&ho_bf, &options, &msh, &stats, &ptree,
         &Ninter, &Nallrows, &Nallcols, &Nalldat_loc,
         allrows, allcols, reinterpret_cast<_Complex float*>(alldat_loc),
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<float>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }
    template<> void bpack_extract_elements_ho<_Complex double>
    (F2Cptr& ho_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, _Complex double* alldat_loc,
     int* rowidx, int* colidx, int* pgidx, int Npmap, int* pmaps) {
      z_c_bpack_extractelement
        (&ho_bf, &options, &msh, &stats, &ptree,
         &Ninter, &Nallrows, &Nallcols, &Nalldat_loc,
         allrows, allcols, reinterpret_cast<_Complex double*>(alldat_loc),
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<double>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }

    template<> void bpack_extract_elements_bf<float>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, float* alldat_loc, int* rowidx, int* colidx,
     int* pgidx, int Npmap, int* pmaps) {
      s_c_bf_extractelement
        (&lr_bf, &options, &msh, &stats, &ptree,
         &Ninter, &Nallrows, &Nallcols, &Nalldat_loc,
         allrows, allcols, alldat_loc,
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<float>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }
    template<> void bpack_extract_elements_bf<double>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, double* alldat_loc, int* rowidx, int* colidx,
     int* pgidx, int Npmap, int* pmaps) {
      d_c_bf_extractelement
        (&lr_bf, &options, &msh, &stats, &ptree,
         &Ninter, &Nallrows, &Nallcols, &Nalldat_loc,
         allrows, allcols, alldat_loc,
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<double>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }
    template<> void bpack_extract_elements_bf<_Complex float>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, _Complex float* alldat_loc,
     int* rowidx, int* colidx, int* pgidx, int Npmap, int* pmaps) {
      c_c_bf_extractelement
        (&lr_bf, &options, &msh, &stats, &ptree,
         &Ninter, &Nallrows, &Nallcols, &Nalldat_loc,
         allrows, allcols, reinterpret_cast<_Complex float*>(alldat_loc),
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<_Complex float>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }
    template<> void bpack_extract_elements_bf<_Complex double>
    (F2Cptr& lr_bf, F2Cptr& options, F2Cptr& msh, F2Cptr& stats,
     F2Cptr& ptree, int Ninter, int Nallrows, int Nallcols, std::int64_t Nalldat_loc,
     int* allrows, int* allcols, _Complex double* alldat_loc,
     int* rowidx, int* colidx, int* pgidx, int Npmap, int* pmaps) {
      z_c_bf_extractelement
        (&lr_bf, &options, &msh, &stats, &ptree,
         &Ninter, &Nallrows, &Nallcols, &Nalldat_loc,
         allrows, allcols, reinterpret_cast<_Complex double*>(alldat_loc),
         rowidx, colidx, pgidx, &Npmap, pmaps);
#if defined(STRUMPACK_COUNT_FLOPS)
      long long int f = bpack_get_stat<_Complex double>(stats, "Flop_C_Extract");
      STRUMPACK_FLOPS(f);
      STRUMPACK_EXTRACTION_FLOPS(f);
#endif
    }

    void bpack_deletestats(F2Cptr& stats) { d_c_bpack_deletestats(&stats); }
    void bpack_deleteproctree(F2Cptr& ptree) { d_c_bpack_deleteproctree(&ptree); }
    void bpack_deletemesh(F2Cptr& mesh) { d_c_bpack_deletemesh(&mesh); }

    void bpack_deletekernelquant(F2Cptr& kerquant) { d_c_bpack_deletekernelquant(&kerquant); }

    void bpack_delete_ho(F2Cptr& ho_bf) { d_c_bpack_delete(&ho_bf); }

    void bpack_delete_bf(F2Cptr& lr_bf) { d_c_bf_deletebf(&lr_bf); }

    void bpack_deleteoptions(F2Cptr& option) { d_c_bpack_deleteoption(&option); }

    template<> void bpack_mult_ho<float>
    (char op, const float* X, float* Y, int Xlrows, int Ylrows, int cols,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      s_c_bpack_mult(&op, X, Y, &Xlrows, &Ylrows,
                     &cols, &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_mult_ho<double>
    (char op, const double* X, double* Y, int Xlrows, int Ylrows, int cols,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      d_c_bpack_mult(&op, X, Y, &Xlrows, &Ylrows,
                     &cols, &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_mult_ho<_Complex float>
    (char op, const _Complex float* X, _Complex float* Y,
     int Xlrows, int Ylrows, int cols, F2Cptr ho_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      c_c_bpack_mult(&op, reinterpret_cast<const _Complex float*>(X),
                     reinterpret_cast<_Complex float*>(Y), &Xlrows, &Ylrows,
                     &cols, &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_mult_ho<_Complex double>
    (char op, const _Complex double* X, _Complex double* Y,
     int Xlrows, int Ylrows, int cols, F2Cptr ho_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      z_c_bpack_mult(&op, reinterpret_cast<const _Complex double*>(X),
                     reinterpret_cast<_Complex double*>(Y), &Xlrows, &Ylrows,
                     &cols, &ho_bf, &options, &stats, &ptree);
    }

    template<> void bpack_mult_bf<float>
    (char op, const float* X, float* Y, int Xlrows, int Ylrows, int cols,
     F2Cptr lr_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      s_c_bf_mult(&op, X, Y, &Xlrows, &Ylrows, &cols,
                  &lr_bf, &options, &stats, &ptree);
    }
    template<> void bpack_mult_bf<double>
    (char op, const double* X, double* Y, int Xlrows, int Ylrows, int cols,
     F2Cptr lr_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      d_c_bf_mult(&op, X, Y, &Xlrows, &Ylrows, &cols,
                  &lr_bf, &options, &stats, &ptree);
    }
    template<> void bpack_mult_bf<_Complex float>
    (char op, const _Complex float* X, _Complex float* Y,
     int Xlrows, int Ylrows, int cols, F2Cptr lr_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      c_c_bf_mult(&op, reinterpret_cast<const _Complex float*>(X),
                  reinterpret_cast<_Complex float*>(Y),
                  &Xlrows, &Ylrows, &cols,
                  &lr_bf, &options, &stats, &ptree);
    }
    template<> void bpack_mult_bf<_Complex double>
    (char op, const _Complex double* X, _Complex double* Y,
     int Xlrows, int Ylrows, int cols, F2Cptr lr_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      z_c_bf_mult(&op, reinterpret_cast<const _Complex double*>(X),
                  reinterpret_cast<_Complex double*>(Y),
                  &Xlrows, &Ylrows, &cols,
                  &lr_bf, &options, &stats, &ptree);
    }


    template<> void bpack_factor_ho<float>
    (F2Cptr& ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree, F2Cptr msh) {
      s_c_bpack_factor(&ho_bf, &options, &stats, &ptree, &msh);
    }
    template<> void bpack_factor_ho<double>
    (F2Cptr& ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree, F2Cptr msh) {
      d_c_bpack_factor(&ho_bf, &options, &stats, &ptree, &msh);
    }
    template<> void bpack_factor_ho<_Complex float>
    (F2Cptr& ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree, F2Cptr msh) {
      c_c_bpack_factor(&ho_bf, &options, &stats, &ptree, &msh);
    }    
    template<> void bpack_factor_ho<_Complex double>
    (F2Cptr& ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree, F2Cptr msh) {
      z_c_bpack_factor(&ho_bf, &options, &stats, &ptree, &msh);
    }


    template<> void bpack_solve_ho<float>
    (float* X, const float* B, int lrows, int rhs,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      s_c_bpack_solve(X, const_cast<float*>(B), &lrows, &rhs,
                      &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_solve_ho<double>
    (double* X, const double* B, int lrows, int rhs,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      d_c_bpack_solve(X, const_cast<double*>(B), &lrows, &rhs,
                      &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_solve_ho<_Complex float>
    (_Complex float* X, const _Complex float* B,
     int lrows, int rhs, F2Cptr ho_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      c_c_bpack_solve
        (reinterpret_cast<_Complex float*>(X),
         reinterpret_cast<_Complex float*>
         (const_cast<_Complex float*>(B)), &lrows, &rhs,
         &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_solve_ho<_Complex double>
    (_Complex double* X, const _Complex double* B,
     int lrows, int rhs, F2Cptr ho_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      z_c_bpack_solve
        (reinterpret_cast<_Complex double*>(X),
         reinterpret_cast<_Complex double*>
         (const_cast<_Complex double*>(B)), &lrows, &rhs,
         &ho_bf, &options, &stats, &ptree);
    }

    template<> void bpack_inv_mult_ho<float>
    (char op, const float* B, float* X, int Xlrows, int Blrows, int rhs,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      s_c_bpack_inv_mult
        (&op, B, X, &Xlrows, &Blrows, &rhs, &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_inv_mult_ho<double>
    (char op, const double* B, double* X, int Xlrows, int Blrows, int rhs,
     F2Cptr ho_bf, F2Cptr options, F2Cptr& stats, F2Cptr ptree) {
      d_c_bpack_inv_mult
        (&op, B, X, &Xlrows, &Blrows, &rhs, &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_inv_mult_ho<_Complex float>
    (char op, const _Complex float* B, _Complex float* X,
     int Xlrows, int Blrows, int rhs, F2Cptr ho_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      c_c_bpack_inv_mult
        (&op, reinterpret_cast<const _Complex float*>(B),
         reinterpret_cast<_Complex float*>(X),
         &Xlrows, &Blrows, &rhs, &ho_bf, &options, &stats, &ptree);
    }
    template<> void bpack_inv_mult_ho<_Complex double>
    (char op, const _Complex double* B, _Complex double* X,
     int Xlrows, int Blrows, int rhs, F2Cptr ho_bf, F2Cptr options,
     F2Cptr& stats, F2Cptr ptree) {
      z_c_bpack_inv_mult
        (&op, reinterpret_cast<const _Complex double*>(B),
         reinterpret_cast<_Complex double*>(X),
         &Xlrows, &Blrows, &rhs, &ho_bf, &options, &stats, &ptree);
    }

    int bpack_treeindex_merged2child(int idx_merge) {
      int idx_child;
      d_c_bpack_treeindex_merged2child(&idx_merge, &idx_child);
      return idx_child;
    }

    int bpack_new2old(F2Cptr msh, int i_new_loc){
      int i_old;
      d_c_bpack_new2old(&msh,&i_new_loc,&i_old);
      return i_old;
    }

    void bpack_getversionnumber(int& v_major, int& v_minor, int& v_bugfix){
      d_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
    }





    // The command line parser for the example related parameters
    void bpack_set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0) {
        
        struct OptionHelp {
            const char* name;
            const char* description;
        };

        static const OptionHelp option_help_table[] = {
          {"nmin_leaf",       "leafsize in the hierarchical partitioning"},
          {"tol_comp",        "relative tolerance for matrix construction"},
          {"tol_rand",        "relative tolerance for matrix inversion"},
          {"tol_Rdetect",     "relative tolerance for rank detection during matrix inversion"},
          {"tol_itersol",     "convergence tolerance for TFQMR iterative solver if precon=2 or 3"},
          {"n_iter",          "maximum iteration count for TFQMR"},
          {"level_check",     "the level in the hierarchical partitioning where the randomized construction algorithm is tested, set to 10000 by default (no checking)"},
          {"precon",          "the use mode of butterflypack: 1: as a direct solver 2: as an iterative solver (compress the matrix and pass it to TFQMR without preconditioner), 3: as a preconditioned iterative solver (compress the matrix and invert the matrix and pass them to TFQMR, using approximate matrix inverse as a preconditioner)"},
          {"xyzsort",         "the hierarchical partitioning algorithm: 0: no permutation 1: permutation based on KD-tree 2: permutation based on cobble-like partitioning"},
          {"lrlevel",         "the level in the hierarchical partitioning (top-down numbered) above which butterfly is used and below which low-rank is used"},
          {"errfillfull",     "errfillfull: a slow (n^2), thorough error checking is performed after the compression of each block"},
          {"baca_batch",      "block size in batched ACA when reclr_leaf=4 or 5"},
          {"reclr_leaf",      "low-rank compression algorithms 1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton 7: ACA with naive parallelization"},
          {"nogeo",           "whether there is geometry information provided 1: is no geometry (xyzsort can not be 1 or 2), 0: there is geometry"},
          {"less_adapt",      "1: improved randomized butterfly construction, default to 1"},
          {"errsol",          "1: generating an artificial true solution vector, compute the RHS with compressed matrix, solve the system, and compare with true solution vector"},
          {"lr_blk_num",      "sqrt of #of subblocks in H-BACA, default to 1"},
          {"rank0",           "initial rank guess in the randomized butterfly algorithm, default to 32"},
          {"rankrate",        "increasing ratio of the rank guess in each iteration, default to 2"},
          {"itermax",         "maximum number of iterations in the randomized butterfly algorithm, default to 10"},
          {"powiter",         "order of power iteration in the randomized low-rank construction"},
          {"ilu",             "whether symmetric gauss-seidel is used when format=2"},
          {"nbundle",         "multiply nbundle sets of vectors together in randomized butterfly algorithm for better flop performance, default to 1"},
          {"near_para",       "admissibility parameter when format=2/3/4/5, strong admissibility typically requires near_para>2.0"},
          {"format",          "the hierarchical matrix format: 1: HODLR/HODBF 2: H matrix 3: HSSBF/SHNBF 4: HSSBF_MD/SHNBF_MD 5: block-LR/BF"},
          {"verbosity",       "verbosity for the printing (-1, 0, 1, 2), -1 suppresses everything, 2 prints most details"},
          {"rmax",            "preestimate of the maximum rank for allocating buffers, default to 1000"},
          {"sample_para",     "oversampling factor in the nlogn entry-evaluation-based butterfly algorithm, default to 2"},
          {"pat_comp",        "pattern of entry-evaluation-based butterfly compression: 1 from right to left, 2 from left to right, 3 from outer to inner"},
          {"knn",             "nearest neighbouring points used in improved BACA and entry-evaluation-based butterfly compression"},
          {"knn_near_para",   "admissibility parameter for guiding the nearest neighbouring points search"},
          {"forwardN15flag",  "whether to use N15 or NlogN algorithm for entry-evaluation-based matrix butterfly compression"},
          {"sample_para_outer","oversampling factor for the outtermost factor matrices in the nlogn entry-evaluation-based butterfly algorithm, default to 2"},
          {"elem_extract",    "0: evaluating entries one by one 1: evaluating entries block by block (may requires communication inside the callback function) 2: evaluating entries block by block (no communication allowed inside the callback function)"},
          {"fastsample_tensor","0: uniformly sample each dimension. 1: uniformly sample the rows of the unfolded matrices on top of 0. 2: use translation invariance"},
          {"use_zfp",         "whether to use zfp compression"},
          {"use_qtt",         "whether to use qtt compression"},
          {"hextralevel",         "HMAT: extra levels for top partitioning of the H matrix based on MPI counts. BLR: Maxlevel-hextralevel is the level for defining B-LR/B-BF blocks"},
          {"help",            "print this help message"}
        };        
                

        double opt_d;
        int opt_i;
        std::vector<std::unique_ptr<char[]>> argv_data(argc);
        std::vector<char*> argv(argc);
        for (int i=0; i<argc; i++) {
          argv_data[i].reset(new char[strlen(cargv[i])+1]);
          argv[i] = argv_data[i].get();
          strcpy(argv[i], cargv[i]);
        }
        option long_options[] =
          {
          {"help",             no_argument,       0, 1000},
          {"nmin_leaf",                     required_argument, 0, 1},  
          {"tol_comp",                   required_argument, 0, 2},  
          {"tol_rand",                   required_argument, 0, 3},  
          {"tol_Rdetect",             required_argument, 0, 4},     
          {"tol_itersol",             required_argument, 0, 5},     
          {"n_iter",          required_argument, 0, 6},			 
          {"level_check",         required_argument, 0, 7},		
          {"precon",                  required_argument, 0, 8}, 	
          {"lrlevel",     required_argument, 0, 10},   
          {"errfillfull",       required_argument, 0, 11}, 
          {"baca_batch",      required_argument, 0, 12}, 
          {"reclr_leaf",      required_argument, 0, 13}, 
          {"nogeo",     required_argument, 0, 14}, 
          {"less_adapt",            required_argument, 0, 15}, 
          {"errsol",           required_argument, 0, 16}, 
          {"lr_blk_num",                  required_argument, 0, 17}, 
          {"rank0",  required_argument, 0, 18}, 
          {"rankrate", required_argument, 0, 19}, 
          {"itermax",               required_argument, 0, 20}, 
          {"powiter",  required_argument, 0, 21}, 
          {"ilu", required_argument, 0, 22},	
          {"nbundle",     required_argument, 0, 23}, 
          {"near_para",  required_argument, 0, 24}, 
          {"format",  required_argument, 0, 25}, 
          {"verbosity", required_argument, 0, 26}, 
          {"rmax", required_argument, 0, 27},	
          {"sample_para", required_argument, 0, 28}, 
          {"pat_comp",    required_argument, 0, 29}, 
          {"knn",         required_argument, 0, 30}, 
          {"knn_near_para",         required_argument, 0, 31}, 
          {"forwardN15flag",         required_argument, 0, 32}, 
          {"sample_para_outer",         required_argument, 0, 33}, 
          {"elem_extract",         required_argument, 0, 34}, 
          {"fastsample_tensor",         required_argument, 0, 35}, 
          {"use_zfp",         required_argument, 0, 36}, 
          {"use_qtt",         required_argument, 0, 37},    
          {"hextralevel",         required_argument, 0, 38},    
          {NULL, 0, NULL, 0}
          };
        int c, option_index = 0;
        // bool unrecognized_options = false;
        opterr = optind = 0;
        while ((c = getopt_long_only
                (argc, argv.data(), "",
                long_options, &option_index)) != -1) {
          switch (c) {
          case 1000: {
            std::cout << "Available ButterflyPACK Command-Line Options:\n";
            for (const auto& opt : option_help_table) {
                std::cout << "  --" << std::setw(20) << std::left << opt.name
                          << " : " << opt.description << "\n";
            }
          } break;
          case 1: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "Nmin_leaf", opt_i);
          } break;
          case 2: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "tol_comp", opt_d);
            d_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
            d_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*0.1);
          } break;
          case 3: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
          } break;
          case 4: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d);
          } break;
          case 5: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "tol_itersol", opt_d);
          } break;
          case 6: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "n_iter", opt_i);
          } break;
          case 7: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "level_check", opt_i);
          } break;
          case 8: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "precon", opt_i);
          } break;
          case 9: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "xyzsort", opt_i);
          } break;
          case 10: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "LRlevel", opt_i);
          } break;
          case 11: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "ErrFillFull", opt_i);
          } break;
          case 12: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "BACA_Batch", opt_i);
          } break;
          case 13: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "RecLR_leaf", opt_i);
          } break;
          case 14: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "nogeo", opt_i);
          } break;
          case 15: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "less_adapt", opt_i);
          } break;
          case 16: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "ErrSol", opt_i);
          } break;
          case 17: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "LR_BLK_NUM", opt_i);
          } break;
          case 18: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "rank0", opt_i);
          } break;
          case 19: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "rankrate", opt_d);
          } break;
          case 20: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "itermax", opt_i);
          } break;
          case 21: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "powiter", opt_i);
          } break;
          case 22: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "ILU", opt_i);
          } break;
          case 23: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "Nbundle", opt_i);
          } break;
          case 24: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "near_para", opt_d);
          } break;
          case 25: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "format", opt_i);
          } break;
          case 26: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "verbosity", opt_i);
          } break;
          case 27: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "rmax", opt_i);
          } break;
          case 28: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "sample_para", opt_d);
          } break;
          case 29: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "pat_comp", opt_i);
          } break;
          case 30: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "knn", opt_i);
          } break;
          case 31: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "knn_near_para", opt_d);
          } break;
          case 32: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "forwardN15flag", opt_i);
          } break;
          case 33: {
            std::istringstream iss(optarg);
            iss >> opt_d;
            d_c_bpack_set_D_option(&option0, "sample_para_outer", opt_d);
          } break;
          case 34: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "elem_extract", opt_i);
          } break;
          case 35: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "fastsample_tensor", opt_i);
          } break;
          case 36: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "use_zfp", opt_i);
          } break;
          case 37: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "use_qtt", opt_i);
          } break;          
          case 38: {
            std::istringstream iss(optarg);
            iss >> opt_i;
            d_c_bpack_set_I_option(&option0, "hextralevel", opt_i);
          } break;
          default: break;
          }
        }
      }

    } // end namespace butterflypack
