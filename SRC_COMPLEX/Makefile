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
SRC = cLAPACK_dgeqp3mod.o cLAPACK_zgeqp3mod.o cLAPACK_sgeqp3mod.o cLAPACK_cgeqp3mod.o cSCALAPACK_pdgeqpfmod.o cSCALAPACK_pzgeqpfmod.o cSCALAPACK_psgeqpfmod.o cSCALAPACK_pcgeqpfmod.o cMISC_linkedlist.o\
cBPACK_defs.o cMISC_DenseLA.o cMISC_utilities.o cBplus_utilities.o cBPACK_utilities.o cBplus_pack_unpack_for_MPI.o\
cBPACK_structure.o cBplus_randomized.o cBplus_compress.o cBplus_factor.o cBPACK_solve_mul.o\
cBPACK_constr.o cBPACK_randomized.o cBPACK_factor.o ccBPACK_wrapper.o

all:lib

lib:$ (SRC)
$ (ARCH) $ (ARCHFLAGS) $ (CButterflyPACKLIB) $ (SRC)
$ (RANLIB) $ (CButterflyPACKLIB)

%.o:%.c
$ (CC) $ (CFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.cpp
$ (CXX) $ (CXXFLAGS) $ (INCDEF) - c$ < $ (VERBOSE)

%.o:%.f
$ (FORTRAN) $ (FFLAGS) - DDAT - c$ < $ (VERBOSE)

%.o:%.f90
$ (FORTRAN) $ (FFLAGS) - DDAT$ (INCDEF) - c$ < $ (VERBOSE)

clean:
rm - f*.o*.mod$ (CButterflyPACKLIB)

