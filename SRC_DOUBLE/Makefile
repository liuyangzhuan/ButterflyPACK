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
SRC = dLAPACK_dgeqp3mod.o dLAPACK_zgeqp3mod.o dLAPACK_sgeqp3mod.o dLAPACK_cgeqp3mod.o dSCALAPACK_pdgeqpfmod.o dSCALAPACK_pzgeqpfmod.o dSCALAPACK_psgeqpfmod.o dSCALAPACK_pcgeqpfmod.o dMISC_linkedlist.o\
dBPACK_defs.o dMISC_DenseLA.o dMISC_utilities.o dBplus_utilities.o dBPACK_utilities.o dBplus_pack_unpack_for_MPI.o\
dBPACK_structure.o dBplus_randomized.o dBplus_compress.o dBplus_factor.o dBPACK_solve_mul.o\
dBPACK_constr.o dBPACK_randomized.o dBPACK_factor.o ddBPACK_wrapper.o

all:lib

lib:$ (SRC)
$ (ARCH) $ (ARCHFLAGS) $ (DButterflyPACKLIB) $ (SRC)
$ (RANLIB) $ (DButterflyPACKLIB)

%.o:%.c
$ (CC) $ (CFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.cpp
$ (CXX) $ (CXXFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.f
$ (FORTRAN) $ (FFLAGS) - DDAT - c$ < $ (VERBOSE)

%.o:%.f90
$ (FORTRAN) $ (FFLAGS) - DDAT$ (INCDEF) - c$ < $ (VERBOSE)

clean:
rm - f*.o*.mod$ (DButterflyPACKLIB)

