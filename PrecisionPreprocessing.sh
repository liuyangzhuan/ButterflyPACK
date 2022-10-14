#!/bin/bash
ROOTDIR=$PWD
SRCDIR=$ROOTDIR/SRC
DSRCDIR=$ROOTDIR/SRC_DOUBLE
ZSRCDIR=$ROOTDIR/SRC_DOUBLECOMPLEX
SSRCDIR=$ROOTDIR/SRC_SINGLE
CSRCDIR=$ROOTDIR/SRC_COMPLEX
CONFIG_FILE=ButterflyPACK_config
MACRO_FILE=$SRCDIR/$CONFIG_FILE.fi
TMP_FILE=$PWD/tmp.txt
lb="\b"
rb="\b"

major=$(grep "set(VERSION_MAJOR" ./CMakeLists.txt | sed 's/[^"]*"\([^"]*\)".*/\1/')
minor=$(grep "set(VERSION_MINOR" ./CMakeLists.txt | sed 's/[^"]*"\([^"]*\)".*/\1/')
bug=$(grep "set(VERSION_BugFix" ./CMakeLists.txt | sed 's/[^"]*"\([^"]*\)".*/\1/')
sed -i -e "s/BPACK_MAJOR_VERSION.*/BPACK_MAJOR_VERSION = $major/" ./SRC/BPACK_defs.f90
sed -i -e "s/BPACK_MINOR_VERSION.*/BPACK_MINOR_VERSION = $minor/" ./SRC/BPACK_defs.f90
sed -i -e "s/BPACK_PATCH_VERSION.*/BPACK_PATCH_VERSION = $bug/" ./SRC/BPACK_defs.f90


######## The following takes care of windows to linux conversion
declare -a StringArray=("*.in" "*.sh" "SRC/*.*" "EXAMPLE/*.*" "Makefile" "*/Makefile")
for val in ${StringArray[@]}; do
   # echo $val
#   sed -i "s/\r$/\r" $val
   sed -i -e "s/[[:blank:]]*$//" $val
done


# sed -i "s/\r$//" *.*
# sed -i "s/\r$//" SRC/*.*
# sed -i "s/\r$//" EXAMPLE/*.*
# sed -i "s/[[:blank:]]*$//" *.*
# sed -i "s/[[:blank:]]*$//" SRC/*.*
# sed -i "s/[[:blank:]]*$//" EXAMPLE/*.*

############################################################################
echo "-- generating macro definition header ..."
rm -rf $TMP_FILE
rm -rf $MACRO_FILE
grep -h "end subroutine" --include='*.f90' --include='*.c' --include='*.h' $SRCDIR/* |sed "s/[[:blank:]]*$//" | sed "s/.* \([^ ][^ ]*\) */\1/" > $TMP_FILE
grep -h "end function" --include='*.f90' --include='*.c' --include='*.h' $SRCDIR/* |sed "s/[[:blank:]]*$//" | sed "s/.* \([^ ][^ ]*\) */\1/" >> $TMP_FILE
grep -h "end type" --include='*.f90' --include='*.c' --include='*.h' $SRCDIR/* |sed "s/[[:blank:]]*$//" | sed "s/.* \([^ ][^ ]*\) */\1/" >> $TMP_FILE
grep -h "end module" --include='*.f90' --include='*.c' --include='*.h' $SRCDIR/* |sed "s/[[:blank:]]*$//" | sed "s/.* \([^ ][^ ]*\) */\1/" >> $TMP_FILE
> $MACRO_FILE
echo "#ifdef DAT" >> $MACRO_FILE
echo "#if DAT==0" >> $MACRO_FILE
echo " " >> $MACRO_FILE
echo "#define DT complex(kind=8)" >> $MACRO_FILE
echo "#define DTR real(kind=8)" >> $MACRO_FILE
echo "#define DTC complex(kind=8)" >> $MACRO_FILE
echo "#define MPI_DT MPI_DOUBLE_COMPLEX" >> $MACRO_FILE
echo "#define C_SIZEOF_DT sizeof_complex16" >> $MACRO_FILE
echo "#define CBIND_DT complex(kind=C_DOUBLE_COMPLEX)" >> $MACRO_FILE
echo "#define C_DT _Complex double" >> $MACRO_FILE
echo "#define gemmf77 zgemm" >> $MACRO_FILE
echo "#define copymatf77 zlacpy" >> $MACRO_FILE
echo " " >> $MACRO_FILE
while IFS= read -r line ; do
    echo "#define $line z_$line" >> $MACRO_FILE
done < "$TMP_FILE"
echo " " >> $MACRO_FILE
echo "#elif DAT==1" >> $MACRO_FILE
echo " " >> $MACRO_FILE
echo "#define DT real(kind=8)" >> $MACRO_FILE
echo "#define DTR real(kind=8)" >> $MACRO_FILE
echo "#define DTC complex(kind=8)" >> $MACRO_FILE
echo "#define MPI_DT MPI_DOUBLE_PRECISION" >> $MACRO_FILE
echo "#define C_SIZEOF_DT sizeof_double" >> $MACRO_FILE
echo "#define CBIND_DT real(kind=C_DOUBLE)" >> $MACRO_FILE
echo "#define C_DT double" >> $MACRO_FILE
echo "#define gemmf77 dgemm" >> $MACRO_FILE
echo "#define copymatf77 dlacpy" >> $MACRO_FILE
echo " " >> $MACRO_FILE
while IFS= read -r line ; do
    echo "#define $line d_$line" >> $MACRO_FILE
done < "$TMP_FILE"
echo " " >> $MACRO_FILE
echo "#elif DAT==2" >> $MACRO_FILE
echo " " >> $MACRO_FILE
echo "#define DT complex(kind=4)" >> $MACRO_FILE
echo "#define DTR real(kind=4)" >> $MACRO_FILE
echo "#define DTC complex(kind=4)" >> $MACRO_FILE
echo "#define MPI_DT MPI_COMPLEX" >> $MACRO_FILE
echo "#define C_SIZEOF_DT sizeof_complex" >> $MACRO_FILE
echo "#define CBIND_DT complex(kind=C_FLOAT_COMPLEX)" >> $MACRO_FILE
echo "#define C_DT _Complex float" >> $MACRO_FILE
echo "#define gemmf77 cgemm" >> $MACRO_FILE
echo "#define copymatf77 clacpy" >> $MACRO_FILE
echo " " >> $MACRO_FILE
while IFS= read -r line ; do
    echo "#define $line c_$line" >> $MACRO_FILE
done < "$TMP_FILE"
echo " " >> $MACRO_FILE
echo "#elif DAT==3" >> $MACRO_FILE
echo " " >> $MACRO_FILE
echo "#define DT real(kind=4)" >> $MACRO_FILE
echo "#define DTR real(kind=4)" >> $MACRO_FILE
echo "#define DTC complex(kind=4)" >> $MACRO_FILE
echo "#define MPI_DT MPI_FLOAT" >> $MACRO_FILE
echo "#define C_SIZEOF_DT sizeof_float" >> $MACRO_FILE
echo "#define CBIND_DT real(kind=C_FLOAT)" >> $MACRO_FILE
echo "#define C_DT float" >> $MACRO_FILE
echo "#define gemmf77 sgemm" >> $MACRO_FILE
echo "#define copymatf77 slacpy" >> $MACRO_FILE
echo " " >> $MACRO_FILE
while IFS= read -r line ; do
    echo "#define $line s_$line" >> $MACRO_FILE
done < "$TMP_FILE"
echo "#endif" >> $MACRO_FILE
echo "#endif" >> $MACRO_FILE
rm -rf $TMP_FILE
cp $MACRO_FILE $ROOTDIR/EXAMPLE/.   # keep a precision-independent ButterflyPACK_config.fi under the example folder

###########################################################
# note that module names and *.h headers need to be renamed without macros
echo "-- copy and modify SRC dir ..."
grep -h "end module" --include='*.f90' --include='*.f' $SRCDIR/* |sed "s/[[:blank:]]*$//" | sed "s/.* \([^ ][^ ]*\) */\1/" > $TMP_FILE


###########################################################
echo "-- processing SRC_DOUBLECOMPLEX ..."
rm -rf $ZSRCDIR
cp -r $SRCDIR $ZSRCDIR
cd $ZSRCDIR
{ echo "#define DAT 0 "; cat $CONFIG_FILE.fi; } >z$CONFIG_FILE.fi
cp z$CONFIG_FILE.fi $ROOTDIR/EXAMPLE/.
for file in *; do
	if [ $file != CMakeLists.txt ] && [ $file != $CONFIG_FILE.fi ] && [ $file != z$CONFIG_FILE.fi ] && [ $file != Makefile ];
	then
		eval sed -i -e "s/$CONFIG_FILE.fi/z$CONFIG_FILE.fi/g" $ZSRCDIR/$file
		eval sed -i -e "s/$file/z$file/g" $ZSRCDIR/CMakeLists.txt
		eval sed -i -e "s/$lb$file$rb/z$file/g" $ZSRCDIR/$file
		objfile=${file%.*}.o
		eval sed -i -e "s/$objfile/z$objfile/g" $ZSRCDIR/Makefile
		mv "$file" "z${file}"
		if [ "$1" = "ON" ];
		then
			cpp -w "z${file}" "z${file}_tmp" # run the cpp preprocessor directly as doxygen can get confused with fortran macros
			mv "z${file}_tmp" "z${file}"
		fi
	fi
done
sed -i -e "s/$CONFIG_FILE.fi/$CONFIG_FILE.fi\\nz$CONFIG_FILE.fi/g" $ZSRCDIR/CMakeLists.txt  # still keep ButterflyPACK_config.fi for backward compatibility
sed -i -e "s/butterflypack/zbutterflypack/g" $ZSRCDIR/CMakeLists.txt
sed -i -e "s/ButterflyPACKLIB/ZButterflyPACKLIB/g" $ZSRCDIR/Makefile
sed -i -e "s/-DDAT/-DDAT=0/g" $ZSRCDIR/CMakeLists.txt
sed -i -e "s/-DDAT/-DDAT=0/g" $ZSRCDIR/Makefile

cd ..
while IFS= read -r line; do
    sed -i -e "s/$lb$line$rb/z_$line/g" $ZSRCDIR/*.f90
    sed -i -e "s/$lb$line$rb/z_$line/g" $ZSRCDIR/*.f
    sed -i -e "s/$lb$line$rb/z_$line/g" $ZSRCDIR/*.h
done < "$TMP_FILE"
sed -i -e "s/C_DT/_Complex double /g" $ZSRCDIR/*.h
sed -i -e "s/c_bpack_/z_c_bpack_/g" $ZSRCDIR/*.h
sed -i -e "s/c_bf_/z_c_bf_/g" $ZSRCDIR/*.h
sed -i -e "s/BPACK_WRAP/z_BPACK_WRAP/g" $ZSRCDIR/*.h
sed -i -e "s/c_bpack_/z_c_bpack_/g" $ZSRCDIR/*.f90
sed -i -e "s/c_bf_/z_c_bf_/g" $ZSRCDIR/*.f90
sed -i -e "s/c_magma_offset_1d/z_c_magma_offset_1d/g" $ZSRCDIR/*.f90
sed -i -e "s/c_magma_offset_2d/z_c_magma_offset_2d/g" $ZSRCDIR/*.f90
sed -i -e "s/magmablas_gemm_vbatched/magmablas_zgemm_vbatched/g" $ZSRCDIR/*.f90


###########################################################
echo "-- processing SRC_DOUBLE ..."
rm -rf $DSRCDIR
cp -r $SRCDIR $DSRCDIR
cd $DSRCDIR
{ echo "#define DAT 1 "; cat $CONFIG_FILE.fi; } >d$CONFIG_FILE.fi
cp d$CONFIG_FILE.fi $ROOTDIR/EXAMPLE/.
for file in *; do
	if [ $file != CMakeLists.txt ] && [ $file != $CONFIG_FILE.fi ] && [ $file != d$CONFIG_FILE.fi ] && [ $file != Makefile ];
	then
		eval sed -i -e "s/$CONFIG_FILE/d$CONFIG_FILE/g" $DSRCDIR/$file
		eval sed -i -e "s/$file/d$file/g" $DSRCDIR/CMakeLists.txt
		eval sed -i -e "s/$lb$file$rb/d$file/g" $DSRCDIR/$file
		objfile=${file%.*}.o
		eval sed -i -e "s/$objfile/d$objfile/g" $DSRCDIR/Makefile
		mv "$file" "d${file}"
		if [ "$1" = "ON" ];
		then
			cpp -w "d${file}" "d${file}_tmp" # run the cpp preprocessor directly as doxygen can get confused with fortran macros
			mv "d${file}_tmp" "d${file}"
		fi
	fi
done
sed -i -e "s/$CONFIG_FILE.fi/$CONFIG_FILE.fi\\nd$CONFIG_FILE.fi/g" $DSRCDIR/CMakeLists.txt   # still keep ButterflyPACK_config.fi for backward compatibility
sed -i -e "s/butterflypack/dbutterflypack/g" $DSRCDIR/CMakeLists.txt
sed -i -e "s/ButterflyPACKLIB/DButterflyPACKLIB/g" $DSRCDIR/Makefile
sed -i -e "s/-DDAT/-DDAT=1/g" $DSRCDIR/CMakeLists.txt
sed -i -e "s/-DDAT/-DDAT=1/g" $DSRCDIR/Makefile

cd ..
while IFS= read -r line; do
	sed -i -e "s/$lb$line$rb/d_$line/g" $DSRCDIR/*.f90
	sed -i -e "s/$lb$line$rb/d_$line/g" $DSRCDIR/*.f
	sed -i -e "s/$lb$line$rb/d_$line/g" $DSRCDIR/*.h
done < "$TMP_FILE"
sed -i -e "s/C_DT/double/g" $DSRCDIR/*.h
sed -i -e "s/c_bpack_/d_c_bpack_/g" $DSRCDIR/*.h
sed -i -e "s/c_bf_/d_c_bf_/g" $DSRCDIR/*.h
sed -i -e "s/BPACK_WRAP/d_BPACK_WRAP/g" $DSRCDIR/*.h
sed -i -e "s/c_bpack_/d_c_bpack_/g" $DSRCDIR/*.f90
sed -i -e "s/c_bf_/d_c_bf_/g" $DSRCDIR/*.f90
sed -i -e "s/c_magma_offset_1d/d_c_magma_offset_1d/g" $DSRCDIR/*.f90
sed -i -e "s/c_magma_offset_2d/d_c_magma_offset_2d/g" $DSRCDIR/*.f90
sed -i -e "s/magmablas_gemm_vbatched/magmablas_dgemm_vbatched/g" $DSRCDIR/*.f90


###########################################################
echo "-- processing SRC_COMPLEX ..."
rm -rf $CSRCDIR
cp -r $SRCDIR $CSRCDIR
cd $CSRCDIR
{ echo "#define DAT 2 "; cat $CONFIG_FILE.fi; } >c$CONFIG_FILE.fi
cp c$CONFIG_FILE.fi $ROOTDIR/EXAMPLE/.
for file in *; do
	if [ $file != CMakeLists.txt ] && [ $file != $CONFIG_FILE.fi ] && [ $file != c$CONFIG_FILE.fi ] && [ $file != Makefile ];
	then
		eval sed -i -e "s/$CONFIG_FILE/c$CONFIG_FILE/g" $CSRCDIR/$file
		eval sed -i -e "s/$file/c$file/g" $CSRCDIR/CMakeLists.txt
		eval sed -i -e "s/$lb$file$rb/c$file/g" $CSRCDIR/$file
		objfile=${file%.*}.o
		eval sed -i -e "s/$objfile/c$objfile/g" $CSRCDIR/Makefile
		mv "$file" "c${file}"
		if [ "$1" = "ON" ];
		then
			cpp -w "c${file}" "c${file}_tmp" # run the cpp preprocessor directly as doxygen can get confused with fortran macros
			mv "c${file}_tmp" "c${file}"
		fi
	fi
done
sed -i -e "s/$CONFIG_FILE.fi/$CONFIG_FILE.fi\\nc$CONFIG_FILE.fi/g" $CSRCDIR/CMakeLists.txt  # still keep ButterflyPACK_config.fi for backward compatibility
sed -i -e "s/butterflypack/cbutterflypack/g" $CSRCDIR/CMakeLists.txt
sed -i -e "s/ButterflyPACKLIB/CButterflyPACKLIB/g" $CSRCDIR/Makefile
sed -i -e "s/-DDAT/-DDAT=2/g" $CSRCDIR/CMakeLists.txt
sed -i -e "s/-DDAT/-DDAT=2/g" $CSRCDIR/Makefile
cd ..
while IFS= read -r line; do
    sed -i -e "s/$lb$line$rb/c_$line/g" $CSRCDIR/*.f90
    sed -i -e "s/$lb$line$rb/c_$line/g" $CSRCDIR/*.f
    sed -i -e "s/$lb$line$rb/c_$line/g" $CSRCDIR/*.h
done < "$TMP_FILE"
sed -i -e "s/C_DT/_Complex float /g" $CSRCDIR/*.h
sed -i -e "s/c_bpack_/c_c_bpack_/g" $CSRCDIR/*.h
sed -i -e "s/c_bf_/c_c_bf_/g" $CSRCDIR/*.h
sed -i -e "s/BPACK_WRAP/c_BPACK_WRAP/g" $CSRCDIR/*.h
sed -i -e "s/c_bpack_/c_c_bpack_/g" $CSRCDIR/*.f90
sed -i -e "s/c_bf_/c_c_bf_/g" $CSRCDIR/*.f90
sed -i -e "s/c_magma_offset_1d/c_c_magma_offset_1d/g" $CSRCDIR/*.f90
sed -i -e "s/c_magma_offset_2d/c_c_magma_offset_2d/g" $CSRCDIR/*.f90
sed -i -e "s/magmablas_gemm_vbatched/magmablas_cgemm_vbatched/g" $CSRCDIR/*.f90


###########################################################
echo "-- processing SRC_SINGLE ..."
rm -rf $SSRCDIR
cp -r $SRCDIR $SSRCDIR
cd $SSRCDIR
{ echo "#define DAT 3 "; cat $CONFIG_FILE.fi; } >s$CONFIG_FILE.fi
cp s$CONFIG_FILE.fi $ROOTDIR/EXAMPLE/.
for file in *; do
	if [ $file != CMakeLists.txt ] && [ $file != $CONFIG_FILE.fi ] && [ $file != s$CONFIG_FILE.fi ] && [ $file != Makefile ];
	then
		eval sed -i -e "s/$CONFIG_FILE/s$CONFIG_FILE/g" $SSRCDIR/$file
		eval sed -i -e "s/$file/s$file/g" $SSRCDIR/CMakeLists.txt
		eval sed -i -e "s/$lb$file$rb/s$file/g" $SSRCDIR/$file
		objfile=${file%.*}.o
		eval sed -i -e "s/$objfile/s$objfile/g" $SSRCDIR/Makefile
		mv "$file" "s${file}"
		if [ "$1" = "ON" ];
		then
			cpp -w "s${file}" "s${file}_tmp" # run the cpp preprocessor directly as doxygen can get confused with fortran macros
			mv "s${file}_tmp" "s${file}"
		fi
	fi
done
sed -i -e "s/$CONFIG_FILE.fi/$CONFIG_FILE.fi\\ns$CONFIG_FILE.fi/g" $SSRCDIR/CMakeLists.txt  # still keep ButterflyPACK_config.fi for backward compatibility
sed -i -e "s/butterflypack/sbutterflypack/g" $SSRCDIR/CMakeLists.txt
sed -i -e "s/ButterflyPACKLIB/SButterflyPACKLIB/g" $SSRCDIR/Makefile
sed -i -e "s/-DDAT/-DDAT=3/g" $SSRCDIR/CMakeLists.txt
sed -i -e "s/-DDAT/-DDAT=3/g" $SSRCDIR/Makefile
cd ..
while IFS= read -r line; do
	sed -i -e "s/$lb$line$rb/s_$line/g" $SSRCDIR/*.f90
	sed -i -e "s/$lb$line$rb/s_$line/g" $SSRCDIR/*.f
	sed -i -e "s/$lb$line$rb/s_$line/g" $SSRCDIR/*.h
done < "$TMP_FILE"
sed -i -e "s/C_DT/float/g" $SSRCDIR/*.h
sed -i -e "s/c_bpack_/s_c_bpack_/g" $SSRCDIR/*.h
sed -i -e "s/c_bf_/s_c_bf_/g" $SSRCDIR/*.h
sed -i -e "s/BPACK_WRAP/s_BPACK_WRAP/g" $SSRCDIR/*.h
sed -i -e "s/c_bpack_/s_c_bpack_/g" $SSRCDIR/*.f90
sed -i -e "s/c_bf_/s_c_bf_/g" $SSRCDIR/*.f90
sed -i -e "s/c_magma_offset_1d/s_c_magma_offset_1d/g" $SSRCDIR/*.f90
sed -i -e "s/c_magma_offset_2d/s_c_magma_offset_2d/g" $SSRCDIR/*.f90
sed -i -e "s/magmablas_gemm_vbatched/magmablas_sgemm_vbatched/g" $SSRCDIR/*.f90




cd $ROOTDIR
rm -rf $TMP_FILE
rm -rf $MACRO_FILE
