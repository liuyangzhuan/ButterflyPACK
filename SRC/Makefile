#######################################################################
#
#  This makefile creates a library for ButterflyPACK.
#  The library is called
#       libbutterflypack.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################
include../make.inc

#
# Routines for ButterflyPACK
#
SRC = LAPACK_dgeqp3mod.o LAPACK_zgeqp3mod.o LAPACK_sgeqp3mod.o LAPACK_cgeqp3mod.o SCALAPACK_pdgeqpfmod.o SCALAPACK_pzgeqpfmod.o SCALAPACK_psgeqpfmod.o SCALAPACK_pcgeqpfmod.o MISC_linkedlist.o\
BPACK_defs.o MISC_DenseLA.o MISC_utilities.o Bplus_utilities.o BPACK_utilities.o Bplus_pack_unpack_for_MPI.o\
BPACK_structure.o Bplus_randomized.o Bplus_compress.o Bplus_factor.o BPACK_solve_mul.o\
BPACK_constr.o BPACK_randomized.o BPACK_factor.o BPACK_wrapper.o

all:lib

lib:$ (SRC)
$ (ARCH) $ (ARCHFLAGS) $ (ButterflyPACKLIB) $ (SRC)
$ (RANLIB) $ (ButterflyPACKLIB)

%.o:%.c
$ (CC) $ (CFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.cpp
$ (CXX) $ (CXXFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.f
$ (FORTRAN) $ (FFLAGS) - DDAT - c$ < $ (VERBOSE)

%.o:%.f90
$ (FORTRAN) $ (FFLAGS) - DDAT$ (INCDEF) - c$ < $ (VERBOSE)

clean:
rm - f*.o*.mod$ (ButterflyPACKLIB)

