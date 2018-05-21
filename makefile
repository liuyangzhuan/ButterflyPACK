# Rules:
.PHONY: all checkdirs clean mvp randh em3d em2d ctest

Compiler=Intel#GNU#
MPI=T
TargetDir = obj
Platform =NERSC#CAC#NERSC#Laptop
GCC=gcc
CFLAGS=-O3 -m64 -openmp


# Flags
ifeq ($(Compiler),Intel)
MODFLAGS = -module $(TargetDir)

ifeq ($(Platform),NERSC)

LIB_MKL = -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64 \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64
INCLUDE_MKL = -I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include/intel64/lp64 -I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include 
F90 = ftn 
CPPC = CC 
LinkFlagC = -dynamic -cxxlib -lifcore  $(LIB_MKL) 
LinkFlagF = -dynamic 
endif

ifeq ($(Platform),Laptop)
LIB_MKL = -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64 \
	-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64
INCLUDE_MKL = -I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include/intel64/lp64 \

F90 = mpiifort
CPPC = icpc 
LinkFlagC = -Bdynamic -cxxlib -lifcore $(LIB_MKL) 
LinkFlagF = -Bdynamic 
endif

ifeq ($(Platform),CAC)
LIB_MKL = -L$(MKL_LIB) \
           -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
INCLUDE_MKL = -I$(MKL_INCLUDE)/intel64/lp64 -I$(MKL_INCLUDE)
F90 = ifort
CPPC = icpc 
LinkFlagC = -Bdynamic -cxxlib -lifcore $(LIB_MKL) 
LinkFlagF = -Bdynamic 
endif

					  
F90FLAGS = -nologo -fpe0 -traceback -cpp -DPRNTlevel=1 -debug full -O0 -g -check bounds -qopenmp -parallel -lpthread $(INCLUDE_MKL) -D$(Compiler)
#F90FLAGS = -O3 -cpp -DPRNTlevel=2 -no-prec-div -axAVX,SSE4.2 -msse2 -align records -parallel -qopenmp -lpthread $(INCLUDE_MKL) -D$(Compiler)  
#CFLAGS= -O0 -g -std=c++11 -qopenmp -debug parallel -traceback 
CFLAGS=-std=c++11 -O3 -qopenmp -qopt-matmul
endif


ifeq ($(Compiler),GNU)
MODFLAGS = -J$(TargetDir) -I$(TargetDir)


ifeq ($(Platform),Laptop)
LIB_MKL = -L/home/administrator/Desktop/software/MKL_INSTALL/lib/intel64 \
       -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread
INCLUDE_MKL = -I/home/administrator/Desktop/software/MKL_INSTALL/include/intel64/lp64 \

F90 = mpif90
CPPC = mpicxx
LinkFlagC = -Bdynamic -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran $(LIB_MKL) 
LinkFlagF = -Bdynamic
endif


ifeq ($(Platform),CAC)
LIB_MKL = -L$(MKL_LIB) \
           -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread
INCLUDE_MKL = -I$(MKL_INCLUDE)/intel64/lp64 -I$(MKL_INCLUDE)
F90 = mpif90
CPPC = mpicxx
LinkFlagC = -Bdynamic -L/usr/lib/gcc/x86_64-linux-gnu/5 -lgfortran $(LIB_MKL) 
LinkFlagF = -Bdynamic
endif
 
F90FLAGS = -DPRNTlevel=2 -O0 -g -pg -cpp -fbacktrace -ffpe-trap=zero,overflow,underflow -fimplicit-none -fbounds-check -ffree-line-length-none  -ffixed-line-length-none -fopenmp -Wconversion -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 $(INCLUDE_MKL)  
#F90FLAGS = -DPRNTlevel=2 -O3 -cpp -ftracer -funswitch-loops -ftree-vectorize -fimplicit-none -fno-range-check -ffree-line-length-none -ffixed-line-length-none -fopenmp -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 $(INCLUDE_MKL) -D$(Compiler) 
CFLAGS=-O0 -g -std=c++11 -fbounds-check -fopenmp -Wconversion -lpthread
#CFLAGS=-std=c++11 -O3 -fopenmp 

endif


SOURCES = Module_file.o misc.o Utilites.o Utilites_randomized.o analytic_part.o gaussian_point.o element_Vinc_cfie.o element_Zmn_cfie.o \
	 	H_matrix_structuring.o Butterfly_compress_forward.o Butterfly_rightmultiply_inverse_randomized.o Randomized_reconstruction.o Bplus_rightmultiply_inverse_randomized.o \
		Butterfly_exact_randomized.o Butterfly_inverse_diagonal_randomized.o Butterfly_inverse_diagonal_randomized_schur_partitioned.o Bplus_inverse_diagonal_randomized_schur_partitioned.o geo_modeling.o current_node_patch_mapping.o RCS_bistatic.o RCS_monostatic.o HODLR_solve.o EM_calculating.o matrices_filling.o HODLR_randomized.o cascading_factorizing.o 
		
obj_Files = $(filter %.o, $(SOURCES))

OBJECTS = $(obj_Files:%.o= $(TargetDir)/%.o)

OBJECTS_fmvpmain = $(TargetDir)/ButterflyMVP.o
OBJECTS_frandhmain = $(TargetDir)/HODLR_Randconst.o
OBJECTS_fem2dmain = $(TargetDir)/EMCURV_Driver.o
OBJECTS_fem3dmain = $(TargetDir)/EMSURF_Driver.o
OBJECTS_cmain = $(TargetDir)/testH.o

all: mvp randh em3d em2d ctest

mvp: EXECUTABLE = fmvpexe
randh: EXECUTABLE = frandhexe
em2d: EXECUTABLE = em2dexe
em3d: EXECUTABLE = em3dexe
ctest: EXECUTABLE = ctestexe

ctest: checkdirs $(OBJECTS) $(OBJECTS_cmain)
	@echo Linking $@ version...
	$(CPPC) $(CFLAGS) $(OBJECTS) $(OBJECTS_cmain) $(LIB_MPI) -o $(EXECUTABLE) $(LinkFlagC)
mvp: checkdirs $(OBJECTS) $(OBJECTS_fmvpmain)
	@echo Linking $@ version...
	$(F90) $(F90FLAGS) $(OBJECTS) $(OBJECTS_fmvpmain) $(LIB_MKL) $(LIB_MPI) -o $(EXECUTABLE) $(LinkFlagF)
randh: checkdirs $(OBJECTS) $(OBJECTS_frandhmain)
	@echo Linking $@ version...
	$(F90) $(F90FLAGS) $(OBJECTS) $(OBJECTS_frandhmain) $(LIB_MKL) $(LIB_MPI) -o $(EXECUTABLE) $(LinkFlagF)
em2d: checkdirs $(OBJECTS) $(OBJECTS_fem2dmain)
	@echo Linking $@ version...
	$(F90) $(F90FLAGS) $(OBJECTS) $(OBJECTS_fem2dmain) $(LIB_MKL) $(LIB_MPI) -o $(EXECUTABLE) $(LinkFlagF)
em3d: checkdirs $(OBJECTS) $(OBJECTS_fem3dmain)
	@echo Linking $@ version...
	$(F90) $(F90FLAGS) $(OBJECTS) $(OBJECTS_fem3dmain) $(LIB_MKL) $(LIB_MPI) -o $(EXECUTABLE) $(LinkFlagF)

	
	
$(TargetDir)/%.o: ./SRC/%.f90 
	$(F90) $(F90FLAGS) $(MODFLAGS) -c $< -o $@
$(TargetDir)/%.o: ./SRC/%.f
	$(F90) $(F90FLAGS) $(MODFLAGS) -c $< -o $@	
$(TargetDir)/%.o: ./SRC/%.c
	$(GCC) $(CFLAGS) -c $< -o $@
$(TargetDir)/%.o: ./SRC/%.cpp
	$(CPPC) $(CFLAGS) -c $< -o $@
	
$(TargetDir)/%.o: ./EXAMPLE/%.f90 
	$(F90) $(F90FLAGS) $(MODFLAGS) -c $< -o $@
$(TargetDir)/%.o: ./EXAMPLE/%.f
	$(F90) $(F90FLAGS) $(MODFLAGS) -c $< -o $@	
$(TargetDir)/%.o: ./EXAMPLE/%.c
	$(GCC) $(CFLAGS) -c $< -o $@
$(TargetDir)/%.o: ./EXAMPLE/%.cpp
	$(CPPC) $(CFLAGS) -c $< -o $@	
	
clean:
	@echo Removing object files in $(TargetDir)/...
	rm -rf $(TargetDir)
	rm -rf ctestexe
	rm -rf fmvpexe
	rm -rf frandhexe
	rm -rf em2dexe
	rm -rf em3dexe
checkdirs: $(TargetDir)

$(TargetDir):
	@echo Creating folder $(TargetDir)/...
	@mkdir -p $@
