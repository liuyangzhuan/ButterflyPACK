######################################################################
#
# CMakeLists.txt for ButterflyPACK
#
######################################################################

# Required version
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

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
set(SUPERLU_VERSION "${PROJECT_VERSION}")
set(SUPERLU_REV "${PROJECT_REV}")


if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

# # setup default options
# set(enable_complex_DEFAULT OFF)


# if (enable_complex)
  # message("-- ButterflyPACK will be built for double-complex data type only.")
  # set(CMAKE_Fortran_FLAGS "-DDAT=0 ${CMAKE_Fortran_FLAGS}")
# else()
  # message("-- ButterflyPACK will be built for double data type only.")
  # set(CMAKE_Fortran_FLAGS "-DDAT=1 ${CMAKE_Fortran_FLAGS}")
# endif()


# # setup required compiler defines and options.
# ## get_directory_property( DirDefs COMPILE_DEFINITIONS )
# # set(CMAKE_C_FLAGS "-DDEBUGlevel=0 -DPRNTlevel=0 ${CMAKE_C_FLAGS}")
# if(XSDK_INDEX_SIZE EQUAL 64)
    # message("-- Using 64 bit integer for index size")
    # set(CMAKE_C_FLAGS "-D_LONGINT ${CMAKE_C_FLAGS}")
# endif()	
# set(CMAKE_C_FLAGS_RELEASE "-O3" CACHE STRING "")

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
#
# # #--------------------- BLAS ---------------------
# if(TPL_BLAS_LIBRARIES)
	# set(BLAS_FOUND TRUE)
# else()
	# find_package(BLAS)
	# if(BLAS_FOUND)
	  # set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH
		# "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
	# endif()
# endif()

# if(BLAS_FOUND)
    # message("-- Using TPL_BLAS_LIBRARIES='${TPL_BLAS_LIBRARIES}'")
    # set(BLAS_LIB ${TPL_BLAS_LIBRARIES})
    # # fix up BLAS library name
    # string (REPLACE ";" " " BLAS_LIB_STR "${BLAS_LIB}")
    # set(BLAS_LIB_EXPORT ${BLAS_LIB_STR})
 # endif()	

 
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
 

#--------------------- MPI ---------------------
find_package(MPI)
if(MPI_Fortran_FOUND )
    set(CMAKE_Fortran_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}" )
endif()

# #--------------------- ParMETIS ---------------------
# if (enable_parmetislib)   ## want to use parmetis
  # if (NOT TPL_PARMETIS_LIBRARIES)
    # message(FATAL_ERROR "TPL_PARMETIS_LIBRARIES option should be set for PARMETIS support to be enabled.")
  # endif()

  # if (NOT TPL_PARMETIS_INCLUDE_DIRS)
    # message(FATAL_ERROR "TPL_PARMETIS_INCLUDE_DIRS option be set for PARMETIS support to be enabled.")
  # endif()
  # foreach(dir ${TPL_PARMETIS_INCLUDE_DIRS})
    # if (NOT EXISTS ${dir})
      # message(FATAL_ERROR "PARMETIS include directory not found: ${dir}")
    # endif()
    # set(CMAKE_C_FLAGS "-I${dir} ${CMAKE_C_FLAGS}")
  # endforeach()

  # message("-- Enabled support for PARMETIS")
  # set(PARMETIS_FOUND TRUE)

  # set(PARMETIS_LIB ${TPL_PARMETIS_LIBRARIES})
  # # fix up PARMETIS library names
  # string (REPLACE ";" " " PARMETIS_LIB_STR "${PARMETIS_LIB}")
  # set(PARMETIS_LIB_EXPORT ${PARMETIS_LIB_STR})

# else()
  # message("-- Will not link with ParMETIS.")
# endif()




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

add_subdirectory(SRC_DOUBLE)
add_subdirectory(SRC_DOUBLECOMPLEX)
add_subdirectory(EXAMPLE)

if(enable_doc)
  message(FATAL_ERROR "Documentation build requested but not implemented.")
  #implement doxygen
endif()

# file(WRITE "make.defs" "# can be exposed to users" ${CMAKE_C_COMPILER})
# configure_file(${CMAKE_SOURCE_DIR}/make.inc.in ${CMAKE_BINARY_DIR}/make.inc)
configure_file(${ButterflyPACK_SOURCE_DIR}/make.inc.in ${ButterflyPACK_SOURCE_DIR}/make.inc)

# Add pkg-config support
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/butterflypack.pc.in ${CMAKE_CURRENT_BINARY_DIR}/butterflypack.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/butterflypack.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
	
	
	
