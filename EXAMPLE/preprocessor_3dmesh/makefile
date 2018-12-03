FC=ifort
#FFLAGS=-O3 -static -ipo -xT #-C -O0 -traceback -g -fpe0 #
#FFLAGS=-vec_report0 -O3 -no-prec-div -xHost -align records
FFLAGS=-O0 -g 
F90SR= nastran2geo.f90
FF90=$(F90SR:.f90=.o)
.SUFFIXES:      .o .f .c .C .h .f90	
.f90.o: 	
	$(FC) $(FFLAGS) -c $<
all:con_all
con_all: $(FF90)
	$(FC) $(FF90) -o $@
clean: 
	rm -f *.o
	rm -f *.mod
	rm -f con_all
