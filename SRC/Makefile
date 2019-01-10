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
include ../make.inc


#
# Routines for ButterflyPACK
#
SRC = dgeqp3mod.o zgeqp3mod.o SCALAPACK_pdgeqpfmod.o SCALAPACK_pzgeqpfmod.o \
	BPACK_defs.o DenseLA.o misc.o Bplus_utilities.o BPACK_utilities.o blocks_pack_unpack_for_MPI.o \
	BPACK_structure.o Bplus_randomized.o Bplus_compress.o Bplus_factor.o BPACK_solve_mul.o \
	BPACK_constr.o BPACK_randomized.o BPACK_factor.o BPACK_wrapper.o

all:  lib

lib: $(SRC)
	$(ARCH) $(ARCHFLAGS) $(ButterflyPACKLIB) $(SRC) 
	$(RANLIB) $(ButterflyPACKLIB)

%.o: %.c
	$(CC) $(CFLAGS) $(INCDEF) -c $< $(VERBOSE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDEF) -c $< $(VERBOSE)

%.o: %.f
	$(FORTRAN) $(FFLAGS) -DDAT -c $< $(VERBOSE)

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -DDAT $(INCDEF) -c $< $(VERBOSE)	

clean:	
	rm -f *.o *.mod $(ButterflyPACKLIB)
	
