############################################################################
#
#  Program:         ButterflyPACK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   January 9, 2019 version 1.0.0
#
#  Modified:  January 9, 2019 version 1.0.0
#
############################################################################

include make.inc

all: lib install example

install:: lib
	( cd SRC_DOUBLE; $(MAKE) )
	( cd SRC_DOUBLECOMPLEX; $(MAKE) )
	( cp $(DButterflyPACKLIB) build/lib; cp $(ZButterflyPACKLIB) build/lib)
	( cp SRC_DOUBLE/*.h build/include; cp SRC_DOUBLECOMPLEX/*.h build/include)

example: lib
	( cd EXAMPLE; $(MAKE) )

clean: cleanlib cleanex

lib:
	( mkdir -p build; cd build; mkdir -p lib; mkdir -p include; mkdir -p SRC_DOUBLE; mkdir -p SRC_DOUBLECOMPLEX)
	( sed -i 's/^M$//' PrecisionPreprocessing.sh; bash PrecisionPreprocessing.sh)
	( cd SRC_DOUBLE; $(MAKE) )
	( cd SRC_DOUBLECOMPLEX; $(MAKE) )

cleanlib:
	( cd SRC_DOUBLE; $(MAKE) clean )
	( cd SRC_DOUBLECOMPLEX; $(MAKE) clean )
	( rm -rf build/lib; rm -rf build/include; rm -rf build/SRC_DOUBLE; rm -rf build/SRC_DOUBLECOMPLEX)

cleanex:
	( cd EXAMPLE; $(MAKE) clean )
