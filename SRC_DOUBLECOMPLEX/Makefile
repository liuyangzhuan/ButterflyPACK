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
SRC = zLAPACK_dgeqp3mod.o zLAPACK_zgeqp3mod.o zLAPACK_sgeqp3mod.o zLAPACK_cgeqp3mod.o zSCALAPACK_pdgeqpfmod.o zSCALAPACK_pzgeqpfmod.o zSCALAPACK_psgeqpfmod.o zSCALAPACK_pcgeqpfmod.o zMISC_linkedlist.o\
zBPACK_defs.o zMISC_DenseLA.o zMISC_utilities.o zBplus_utilities.o zBPACK_utilities.o zBplus_pack_unpack_for_MPI.o\
zBPACK_structure.o zBplus_randomized.o zBplus_compress.o zBplus_factor.o zBPACK_solve_mul.o\
zBPACK_constr.o zBPACK_randomized.o zBPACK_factor.o zzBPACK_wrapper.o

all:lib

lib:$ (SRC)
$ (ARCH) $ (ARCHFLAGS) $ (ZButterflyPACKLIB) $ (SRC)
$ (RANLIB) $ (ZButterflyPACKLIB)

%.o:%.c
$ (CC) $ (CFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.cpp
$ (CXX) $ (CXXFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.f
$ (FORTRAN) $ (FFLAGS) - DDAT - c$ < $ (VERBOSE)

%.o:%.f90
$ (FORTRAN) $ (FFLAGS) - DDAT$ (INCDEF) - c$ < $ (VERBOSE)

clean:
rm - f*.o*.mod$ (ZButterflyPACKLIB)

