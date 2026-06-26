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
SRC = sLAPACK_dgeqp3mod.o sLAPACK_zgeqp3mod.o sLAPACK_sgeqp3mod.o sLAPACK_cgeqp3mod.o sSCALAPACK_pdgeqpfmod.o sSCALAPACK_pzgeqpfmod.o sSCALAPACK_psgeqpfmod.o sSCALAPACK_pcgeqpfmod.o sMISC_linkedlist.o\
sBPACK_defs.o sMISC_DenseLA.o sMISC_utilities.o sBplus_utilities.o sBPACK_utilities.o sBplus_pack_unpack_for_MPI.o\
sBPACK_structure.o sBplus_randomized.o sBplus_compress.o sBplus_factor.o sBPACK_solve_mul.o\
sBPACK_constr.o sBPACK_randomized.o sBPACK_factor.o ssBPACK_wrapper.o

all:lib

lib:$ (SRC)
$ (ARCH) $ (ARCHFLAGS) $ (SButterflyPACKLIB) $ (SRC)
$ (RANLIB) $ (SButterflyPACKLIB)

%.o:%.c
$ (CC) $ (CFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.cpp
$ (CXX) $ (CXXFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.f
$ (FORTRAN) $ (FFLAGS) - DDAT - c$ < $ (VERBOSE)

%.o:%.f90
$ (FORTRAN) $ (FFLAGS) - DDAT$ (INCDEF) - c$ < $ (VERBOSE)

clean:
rm - f*.o*.mod$ (SButterflyPACKLIB)

