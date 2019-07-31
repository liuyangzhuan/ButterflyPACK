######################################################################
#
# CMakeLists.txt for ButterflyPACK
#
######################################################################

# Required version
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# SET(CMAKE_INSTALL_PREFIX < install_path >)

# Project version numbers
project(ButterflyPACK C CXX Fortran)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_BugFix "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# IDEAS: xSDK standards module
MESSAGE("\nProcess XSDK defaults ...")
# SET(USE_XSDK_DEFAULTS_DEFAULT TRUE) # Set to false if desired
INCLUDE("cmake/XSDKDefaults.cmake")
INCLUDE(CTest)

######################################################################
#
# Load all macros:
#
FILE(GLOB _macro_files "cmake/Macros/*.cmake")
MESSAGE(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_external_macros.cmake")
FOREACH(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

######################################################################
#
# Usual initialization stuff
#
######################################################################
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)    ## ????
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

#---- For shared library

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#----

if (BUILD_SHARED_LIBS)
  message("-- ButterflyPACK will be built as a shared library.")
else()
  message("-- ButterflyPACK will be built as a static library.")
endif()

enable_language (C)
enable_language (CXX)
set(CMAKE_CXX_STANDARD 11)
if (XSDK_ENABLE_Fortran)
  enable_language (Fortran)
  set(NOFORTRAN FALSE)
endif()
set(ButterflyPACK_VERSION "${PROJECT_VERSION}")
set(ButterflyPACK_REV "${PROJECT_REV}")


if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()


######################################################################
#
# Add compiler-specific compiler flags
#
######################################################################

set(CMAKE_Fortran_FLAGS "-cpp ${CMAKE_Fortran_FLAGS}")
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-ffree-line-length-none" COMPILER_GNU)
if (COMPILER_GNU)
set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -ffixed-line-length-none -fimplicit-none -lpthread ${CMAKE_Fortran_FLAGS}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
set(CMAKE_Fortran_FLAGS "-fbacktrace -fbounds-check -Wconversion ${CMAKE_Fortran_FLAGS}")
endif()
endif()
check_fortran_compiler_flag("-no-prec-div" COMPILER_Intel)
if (COMPILER_Intel)
set(CMAKE_Fortran_FLAGS "-DIntel ${CMAKE_Fortran_FLAGS}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
set(CMAKE_Fortran_FLAGS "-traceback -debug full -check bounds ${CMAKE_Fortran_FLAGS}")
endif()
endif()


set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-qopt-matmul" qoptmatmul)
if (qoptmatmul)
set(CMAKE_CXX_FLAGS  "-qopt-matmul ${CMAKE_CXX_FLAGS}")
endif()



######################################################################
#
# Find packages
#
######################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#--------------------- OpenMP ---------------------
find_package(OpenMP)
## include(FindOpenMP)  # Strumpack uses this
if(OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS}")
# On edison, OpenMP_EXE_LINKER_FLAGS is empty
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")
 message("-- OpenMP_EXE_LINKER_FLAGS='${OpenMP_EXE_LINKER_FLAGS}'")
 message("-- CMAKE_EXE_LINKER_FLAGS='${CMAKE_EXE_LINKER_FLAGS}'")
endif()



 # #--------------------- BLAS ---------------------
if(TPL_BLAS_LIBRARIES)
	set(BLAS_FOUND TRUE)
else()
	find_package(BLAS)
	if(BLAS_FOUND)
	  set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH
		"Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
	endif()
endif()

if(BLAS_FOUND)
    message("-- Using TPL_BLAS_LIBRARIES='${TPL_BLAS_LIBRARIES}'")
    set(BLAS_LIB ${TPL_BLAS_LIBRARIES})
    # fix up BLAS library name
    string (REPLACE ";" " " BLAS_LIB_STR "${BLAS_LIB}")
    set(BLAS_LIB_EXPORT ${BLAS_LIB_STR})
 endif()


 # #--------------------- LAPACK ---------------------
if(TPL_LAPACK_LIBRARIES)
	set(LAPACK_FOUND TRUE)
else()
	find_package(LAPACK)
	if(LAPACK_FOUND)
	  set(TPL_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}" CACHE FILEPATH
		"Set from FindLAPACK.cmake LAPACK_LIBRARIES." FORCE)
	endif()
endif()

if(LAPACK_FOUND)
    message("-- Using TPL_LAPACK_LIBRARIES='${TPL_LAPACK_LIBRARIES}'")
    set(LAPACK_LIB ${TPL_LAPACK_LIBRARIES})
    # fix up LAPACK library name
    string (REPLACE ";" " " LAPACK_LIB_STR "${LAPACK_LIB}")
    set(LAPACK_LIB_EXPORT ${LAPACK_LIB_STR})
 endif()

  # #--------------------- SCALAPACK ---------------------
if(TPL_SCALAPACK_LIBRARIES)
	set(SCALAPACK_FOUND TRUE)
else()
	find_package(SCALAPACK)
	if(SCALAPACK_FOUND)
	  set(TPL_SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARIES}" CACHE FILEPATH
		"Set from FindSCALAPACK.cmake SCALAPACK_LIBRARIES." FORCE)
	endif()
endif()

if(SCALAPACK_FOUND)
    message("-- Using TPL_SCALAPACK_LIBRARIES='${TPL_SCALAPACK_LIBRARIES}'")
    set(SCALAPACK_LIB ${TPL_SCALAPACK_LIBRARIES})
    # fix up LAPACK library name
    string (REPLACE ";" " " SCALAPACK_LIB_STR "${SCALAPACK_LIB}")
    set(SCALAPACK_LIB_EXPORT ${SCALAPACK_LIB_STR})
 endif()


  # #--------------------- ARPACK ---------------------

if(TPL_ARPACK_LIBRARIES)
	set(ARPACK_FOUND TRUE)
else()
	find_package(ArpackNg)
	if(SCALAPACK_FOUND)
	  set(TPL_SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARIES}" CACHE FILEPATH
		"Set from FindSCALAPACK.cmake SCALAPACK_LIBRARIES." FORCE)
	endif()
endif()

if(TPL_ARPACK_LIBRARIES)
	set(CMAKE_Fortran_FLAGS "-DHAVE_ARPACK ${CMAKE_Fortran_FLAGS}")
    message("-- Using TPL_ARPACK_LIBRARIES='${TPL_ARPACK_LIBRARIES}'")
    set(ARPACK_LIB ${TPL_ARPACK_LIBRARIES})
    # fix up LAPACK library name
    string (REPLACE ";" " " ARPACK_LIB_STR "${ARPACK_LIB}")
    set(ARPACK_LIB_EXPORT ${ARPACK_LIB_STR})
endif()



#--------------------- MPI ---------------------
find_package(MPI)
if(MPI_Fortran_FOUND )
    set(CMAKE_Fortran_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}" )
endif()

#---------------------- Additional Fortran linker library ---------
SET(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
# SET_IF_EMPTY(_fortran_libs gfortran m quadmath c)
FOREACH(_lib ${_fortran_libs})
  # FIND_SYSTEM_LIBRARY(${_lib}_LIBRARY NAMES ${_lib})
  set(EXTRA_LIB "-l${_lib} ${EXTRA_LIB}")
ENDFOREACH()



######################################################################
#
# Include directories
#
######################################################################

# include_directories(${CMAKE_SOURCE_DIR}/SRC)
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${TPL_LAPACK95_INCLUDE_DIRS})

######################################################################
#
# Add subdirectories
#
######################################################################

execute_process (
    COMMAND bash PrecisionPreprocessing.sh
	WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
)

add_subdirectory(SRC_DOUBLE)
add_subdirectory(SRC_DOUBLECOMPLEX)
add_subdirectory(EXAMPLE)

if(enable_doc)
  message(FATAL_ERROR "Documentation build requested but not implemented.")
  #implement doxygen
endif()

configure_file(${ButterflyPACK_SOURCE_DIR}/make.inc.in ${ButterflyPACK_SOURCE_DIR}/make.inc)

# Add pkg-config support
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/butterflypack.pc.in ${CMAKE_CURRENT_BINARY_DIR}/butterflypack.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/butterflypack.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

