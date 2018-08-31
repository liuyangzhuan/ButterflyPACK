############################################################################
#
#  Program:         HODLRBF
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   August 19, 2018	version 1.0.0
#
#  Modified:	    
#		    
#
############################################################################
#
#  The name of the libraries to be created/linked to
#
HODLRBFroot	= /home/administrator/Desktop/research/hod-lr-bf/build
DHODLRBFLIB   	= $(HODLRBFroot)/SRC_DOUBLE/libdhodlrbf.a $(HODLRBFroot)/SRC_DOUBLECOMPLEX/libzhodlrbf.a

LIBS		= $(DHODLRBFLIB) /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_gf_lp64.so /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_intel_thread.so /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_core.so /opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/libiomp5.so /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so /opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64/libmkl_scalapack_lp64.so

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = /usr/bin/ranlib

CC           = mpicxx ##/opt/cray/pe/craype/2.5.12/bin/cc
CFLAGS 	     = -O3 -DNDEBUG 
# CFLAGS       += -D_LONGINT   ## 64-bit integer
# CFLAGS       += -D
# CFLAGS       +=  
NOOPTS       = -O0
FORTRAN	     = mpif90 ## /opt/cray/pe/craype/2.5.12/bin/ftn

LOADER       = $(CC)
LOADOPTS     = -Wl,-rpath,-qopenmp  -qopenmp 
