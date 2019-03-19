#######################################################################
#
#  This makefile creates the example programs for ButterflyPACK.
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#
#######################################################################
include ../make.inc

FZEXM1 = EMCURV_Module.o EMCURV_Driver.o  
FZEXM2 = EMSURF_Module.o EMSURF_Driver.o  
FZEXM3 = SMAT_Driver.o
FZEXM4 = EMCURV_Module.o EMCURV_Eigen_Driver.o
FZEXM5 = FrontalDist_Driver.o
FZEXM6 = EMSURF_Module.o EMSURF_Eigen_Driver.o 

FDEXM1 = FULLMAT_Driver.o
FDEXM2 = KERREG_Driver.o
# FDEXM3 = FULLMATKERREG_Driver.o
FDEXM4 = Frontal_Driver.o

CDEXM1 = InterfaceTest.o

all: double doublecomplex

double: ctest full krr frontal frontaldist
# double: ctest full krr fullkrr frontal

doublecomplex: ie2d ie3d smat ie2deigen ie3deigen 
 
ie2d: $(FZEXM1) $(ZButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FZEXM1) $(LIBS) -lm -o $@
	
ie3d: $(FZEXM2) $(ZButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FZEXM2) $(LIBS) -lm -o $@

smat: $(FZEXM3) $(ZButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FZEXM3) $(LIBS) -lm -o $@

ie2deigen: $(FZEXM4) $(ZButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FZEXM4) $(LIBS) -lm -o $@	
	
frontaldist: $(FZEXM5) $(ZButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FZEXM5) $(LIBS) -lm -o $@	
	
ie3deigen: $(FZEXM6) $(ZButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FZEXM6) $(LIBS) -lm -o $@	
	
	
full: $(FDEXM1) $(DButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FDEXM1) $(LIBS) -lm -o $@

krr: $(FDEXM2) $(DButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FDEXM2) $(LIBS) -lm -o $@	
	
# fullkrr: $(FDEXM3) $(DButterflyPACKLIB)
	# $(FLOADER) $(FLOADOPTS) $(FDEXM3) $(LIBS) -lm -o $@		
	
frontal: $(FDEXM4) $(DButterflyPACKLIB)
	$(FLOADER) $(FLOADOPTS) $(FDEXM4) $(LIBS) -lm -o $@		
	
ctest: $(CDEXM1) $(DButterflyPACKLIB)
	$(CXXLOADER) $(CXXLOADOPTS) $(CDEXM1) $(LIBS) -lm $(FIMPLICIT_LINKLIB) -o $@		
	

%.o: %.c
	$(CC) $(CFLAGS) $(INCDEF) -c $< $(VERBOSE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDEF) -c $< $(VERBOSE)

%.o: %.f
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

%.o: %.f90
	$(FORTRAN) $(FFLAGS) $(INCDEF) -c $< $(VERBOSE)		
	

clean:	
	# rm -f *.o *.mod ie2d ie3d smat ctest full krr fullkrr
	rm -f *.o *.mod ie2d ie3d smat ctest full krr frontal ie2deigen ie3deigen frontaldist


