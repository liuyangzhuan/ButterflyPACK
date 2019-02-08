# ButterflyPACK
ButterflyPACK, Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

[![Build Status](https://travis-ci.com/liuyangzhuan/ButterflyPACK.svg?token=xooeQZbwgfe8y48ztwEU&branch=master)](https://travis-ci.com/liuyangzhuan/ButterflyPACK)


## Overview
ButterflyPACK is a mathematical software for rapidly solving large-scale dense linear systems that exhibit off-diagonal rank-deficiency. These systems arise frequently from boundary element methods, or factorization phases in finite-differenece/finite-element methods. ButterflyPACK relies on low-rank or butterfly formats under Hierarchical matrix, HODLR or other hierarhically nested frameworks to compress, factor and solve the linear system in quasi-linear time. The computationally most intensive phase, factorization, is accelerated via randomized linear algebras. The butterfly format, originally inspired by the butterfly data flow in fast Fourier Transform, is a linear algebra tool well-suited for compressing matrices arising from high-frequency wave equations or highly oscillatory integral operators. ButterflyPACK also provides preconditioned TFQMR iterative solvers.

ButterflyPACK is written in Fortran 2003, it also has C++ interfaces. ButterflyPACK supports hybrid MPI/OpenMP programming models. In addition, ButterflyPACK can be readily invoked from the software STRUMPACK for solving dense and sparse linear systems.


## Installation

The installation uses CMake build system. You may need "dos2unix" and "bash" for the build process. The software also requires BLAS, LAPACK and SCALAPACK.

For an installation with GNU compiliers, do:
```
export BLAS_LIB=<Lib of the BLAS installation>
export LAPACK_LIB=<Lib of the LAPACK installation>
export SCALAPACK_LIB=<Lib of the SCALAPACK installation>
sh PrecisionPreprocessing.sh
mkdir build ; cd build;
cmake .. \
	-DCMAKE_Fortran_FLAGS="" \
	-DCMAKE_CXX_FLAGS="" \
	-DTPL_BLAS_LIBRARIES="${BLAS_LIB}" \
	-DTPL_LAPACK_LIBRARIES="${LAPACK_LIB}" \
	-DTPL_SCALAPACK_LIBRARIES="${SCALAPACK_LIB}" \
	-DBUILD_SHARED_LIBS=ON \
	-DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_C_COMPILER=mpicc \
	-DCMAKE_INSTALL_PREFIX=. \
	-DCMAKE_BUILD_TYPE=Release
make
( see example cmake script in example_scripts: run_cmake_build_gnu_ubuntu.sh, run_cmake_build_intel_ubuntu.sh, run_cmake_build_gnu_cori.sh, run_cmake_build_intel_cori.sh)
```

## Current developers
 - Yang Liu - liuyangzhuan@lbl.gov (Lawrence Berkeley National Laboratory)

## Other contributors
 - Wissam Sid-Lakhdar - wissam@lbl.gov (Lawrence Berkeley National Laboratory)
 - Pieter Ghysels - pghysels@lbl.gov (Lawrence Berkeley National Laboratory)
 - Xiaoye S. Li - xsli@lbl.gov (Lawrence Berkeley National Laboratory)
 - Han Guo - hanguo@umich.edu (University of Michigan)
 - Eric Michielssen - emichiel@umich.edu (University of Michigan)
 - Haizhao Yang - matyh@nus.edu.sg (National University of Singapore)

## Reference
to be added ...
