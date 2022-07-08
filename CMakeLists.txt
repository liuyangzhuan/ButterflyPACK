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
set(VERSION_MAJOR "2")
set(VERSION_MINOR "1")
set(VERSION_BugFix "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

# Set up options
option(enable_doc       "Build doxygen documentation" OFF)
option(enable_openmp    "Build with OpenMP support" ON)

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
SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

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

if(BUILD_SHARED_LIBS)
  message("-- ButterflyPACK will be built as a shared library.")
else()
  message("-- ButterflyPACK will be built as a static library.")
endif()

enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
if(XSDK_ENABLE_Fortran)
  enable_language(Fortran)
  set(NOFORTRAN FALSE)
endif()
set(ButterflyPACK_VERSION "${PROJECT_VERSION}")
set(ButterflyPACK_REV "${PROJECT_REV}")

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

######################################################################
#
# Add compiler-specific compiler flags
#
######################################################################
include(CheckFortranCompilerFlag)
include(CheckCXXCompilerFlag)

check_fortran_compiler_flag("-ffree-line-length-none" COMPILER_GNU)
if(COMPILER_GNU)
  set(CMAKE_Fortran_FLAGS "-cpp -fno-range-check -ffree-line-length-none -ffixed-line-length-none -fimplicit-none -lpthread ${CMAKE_Fortran_FLAGS}")
  check_fortran_compiler_flag("-fallow-argument-mismatch" GNU10)
  if(GNU10)
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch ${CMAKE_Fortran_FLAGS}")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "-fno-range-check -fbacktrace -fbounds-check -Wconversion ${CMAKE_Fortran_FLAGS}")
  endif()
endif()

check_fortran_compiler_flag("-no-prec-div" COMPILER_Intel)
if(COMPILER_Intel)
  set(CMAKE_Fortran_FLAGS "-cpp -DIntel ${CMAKE_Fortran_FLAGS}")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "-traceback -debug full -check bounds ${CMAKE_Fortran_FLAGS}")
  endif()
endif()

check_fortran_compiler_flag("-Mextend" COMPILER_PGI)
if(COMPILER_PGI)
  set(CMAKE_Fortran_FLAGS "-cpp -Mextend -DPGI ${CMAKE_Fortran_FLAGS}")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "-O0 -g -Mbounds -Mdepchk -traceback ${CMAKE_Fortran_FLAGS}")
  endif()
endif()

check_fortran_compiler_flag("-N 1023" COMPILER_CRAY)
if(COMPILER_CRAY OR CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS "-e Z -h omp -N 1023 -DCRAY ${CMAKE_Fortran_FLAGS}")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "-O0 -g -e D ${CMAKE_Fortran_FLAGS}")
  endif()
endif()

# Test OPENMP TASKLOOP
if(enable_openmp)
  if(COMPILER_PGI OR COMPILER_CRAY OR CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")  # PGI and Cray taskloop seems buggy
  else()
    check_fortran_source_compiles(
"PROGRAM TEST
integer i
!$omp parallel
!$omp single
!$omp taskloop default(shared) private(i)
do i=1,10
write(*,*)'iter:',i
enddo
!$omp end taskloop
!$omp end single
!$omp end parallel
END PROGRAM TEST" BUTTERFLYPACK_USE_TASKLOOP SRC_EXT F90)
    if(BUTTERFLYPACK_USE_TASKLOOP)
      set(CMAKE_Fortran_FLAGS "-DHAVE_TASKLOOP ${CMAKE_Fortran_FLAGS}")
      message("-- Using OpenMP taskloop")
    endif()
  endif()
endif()

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
check_cxx_compiler_flag("-qopt-matmul" qoptmatmul)
if(qoptmatmul)
  set(CMAKE_CXX_FLAGS  "-qopt-matmul ${CMAKE_CXX_FLAGS}")
endif()

######################################################################
#
# Find packages
#
######################################################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#-------- Additional Fortran linker library --------
#if(APPLE)
#else()
set(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
foreach(_lib ${_fortran_libs})
  # FIND_SYSTEM_LIBRARY(${_lib}_LIBRARY NAMES ${_lib})
  set(EXTRA_LIB "-l${_lib} ${EXTRA_LIB}")
endforeach()
#endif()

#--------------------- OpenMP ---------------------
if(enable_openmp)
  find_package(OpenMP)
  if(OPENMP_FORTRAN_FOUND)
    set(CMAKE_Fortran_FLAGS "-DHAVE_OPENMP ${OpenMP_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS}")
    set(OpenMP_Fortran_FLAGS_EXPORT "${OpenMP_Fortran_FLAGS}")
    set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ${OpenMP_Fortran_FLAGS}")
  endif()
  if(OPENMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "-DHAVE_OPENMP ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(OpenMP_CXX_FLAGS_EXPORT "${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

#--------------------- BLAS ---------------------
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
  string(REPLACE ";" " " BLAS_LIB_STR "${BLAS_LIB}")
  set(BLAS_LIB_EXPORT ${BLAS_LIB_STR})
  string(FIND ${BLAS_LIB_STR} "mkl" BUTTERFLYPACK_USE_MKL)
  if(BUTTERFLYPACK_USE_MKL GREATER -1)
    set(CMAKE_Fortran_FLAGS "-DHAVE_MKL ${CMAKE_Fortran_FLAGS}")
    message("-- Using MKL VSL and batch GEMM")
  endif()
endif()

#--------------------- LAPACK ---------------------
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
  string(REPLACE ";" " " LAPACK_LIB_STR "${LAPACK_LIB}")
  set(LAPACK_LIB_EXPORT ${LAPACK_LIB_STR})
endif()

#--------------------- SCALAPACK ---------------------
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
  # fix up ScaLAPACK library name
  string(REPLACE ";" " " SCALAPACK_LIB_STR "${SCALAPACK_LIB}")
  set(SCALAPACK_LIB_EXPORT ${SCALAPACK_LIB_STR})
endif()

#--------------------- MAGMA ---------------------
if(TPL_MAGMA_LIBRARIES)
  set(MAGMA_FOUND TRUE)
endif()

if(MAGMA_FOUND)
  message("-- Using TPL_MAGMA_LIBRARIES='${TPL_MAGMA_LIBRARIES}'")
  set(MAGMA_LIB ${TPL_MAGMA_LIBRARIES})
  # fix up MAGMA library name
  string(REPLACE ";" " " MAGMA_LIB_STR "${MAGMA_LIB}")
  set(MAGMA_LIB_EXPORT ${MAGMA_LIB_STR})
  set(CMAKE_Fortran_FLAGS "-DHAVE_MAGMA ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "-DHAVE_MAGMA ${CMAKE_C_FLAGS}")
endif()

#--------------------- ARPACK ---------------------
if(TPL_ARPACK_LIBRARIES)
  set(CMAKE_Fortran_FLAGS "-DHAVE_ARPACK ${CMAKE_Fortran_FLAGS}")
  message("-- Using TPL_ARPACK_LIBRARIES='${TPL_ARPACK_LIBRARIES}'")
  set(ARPACK_LIB ${TPL_ARPACK_LIBRARIES})
  # fix up ARPACK library name
  string(REPLACE ";" " " ARPACK_LIB_STR "${ARPACK_LIB}")
  set(ARPACK_LIB_EXPORT ${ARPACK_LIB_STR})
endif()

#--------------------- MPI ---------------------
find_package(MPI)
if(MPI_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
endif()

# Test MPI3
check_cxx_source_compiles(
"#include <mpi.h>
#include <iostream>
int main(int argc, char *argv[]) {
int x=0;
MPI_Request req;

MPI_Init(&argc,&argv);
MPI_Ibcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD,&req);
MPI_Finalize();

return 0;
}" BUTTERFLYPACK_USE_MPI3)

if(BUTTERFLYPACK_USE_MPI3)
  set(CMAKE_Fortran_FLAGS "-DHAVE_MPI3 ${CMAKE_Fortran_FLAGS}")
  message("-- Using MPI3 features")
endif()

# Test Fortran Destructor
check_fortran_source_compiles(
"module MOD
  type quant
    contains
    final :: finalize
  end type quant

  contains

  subroutine finalize( this )
    type(quant), intent(inout) :: this
  end subroutine finalize
end module MOD

PROGRAM TEST
  use MOD
END PROGRAM TEST" BUTTERFLYPACK_USE_Finalizer SRC_EXT F90)

if(BUTTERFLYPACK_USE_Finalizer)
  set(CMAKE_Fortran_FLAGS "-DHAVE_FINAL ${CMAKE_Fortran_FLAGS}")
  message("-- Using Fortran Finalizer features")
endif()

######################################################################
#
# Include directories
#
######################################################################
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${TPL_LAPACK95_INCLUDE_DIRS})

######################################################################
#
# Add subdirectories
#
######################################################################
find_program(SED_TOOL NAMES sed)
if(NOT SED_TOOL)
  message(FATAL_ERROR "Unable to find sed")
endif()

if(APPLE)
  # By default assume BSD sed command is used on macOS
  set(USING_MACOS_SED TRUE)
  set(SED_ARG "\"\"") # mac (BSD) sed

  # Detect version of SED_TOOL version
  message(STATUS "detecting sed version...")
  execute_process(
    COMMAND ${SED_TOOL} --version
    WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    OUTPUT_VARIABLE SED_TOOL_VERSION
    ERROR_QUIET
  )
  message(STATUS "${SED_TOOL_VERSION}")

  if(SED_TOOL_VERSION MATCHES ".*GNU sed.*")
      message(STATUS "using GNU sed...")
      set(SED_ARG "") # GNU sed
      set(USING_MACOS_SED FALSE)
  endif()

  message(STATUS "sed command [${SED_TOOL} -i ${SED_ARG}]")

  execute_process(
    COMMAND cp PrecisionPreprocessing.sh PrecisionPreprocessing_mac.sh
    WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
  )

  if(USING_MACOS_SED)
    execute_process(
      COMMAND ${SED_TOOL} -i ${SED_ARG} -e "s/sed -i/sed -i \"\"/g" PrecisionPreprocessing_mac.sh
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    execute_process(
      COMMAND ${SED_TOOL} -i ${SED_ARG} -e "s/eval//g" PrecisionPreprocessing_mac.sh
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    execute_process(
      COMMAND ${SED_TOOL} -i ${SED_ARG} -e "/\r\"/d" PrecisionPreprocessing_mac.sh
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    execute_process(
      COMMAND ${SED_TOOL} -i ${SED_ARG} -e "s/lb=.*/lb=\"[[:<:]]\"/g" PrecisionPreprocessing_mac.sh
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    execute_process(
      COMMAND ${SED_TOOL} -i ${SED_ARG} -e "s/rb=.*/rb=\"[[:>:]]\"/g" PrecisionPreprocessing_mac.sh
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    execute_process(
      COMMAND bash PrecisionPreprocessing_mac.sh ${enable_doc}
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    set(PREPROC_EXPORT "PrecisionPreprocessing_mac.sh")
  else()
    execute_process(
      COMMAND bash PrecisionPreprocessing.sh ${enable_doc}
      WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
    )
    set(PREPROC_EXPORT "PrecisionPreprocessing.sh")
  endif()
elseif(UNIX)
  execute_process(
    COMMAND bash PrecisionPreprocessing.sh ${enable_doc}
    WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
  )
  set(PREPROC_EXPORT "PrecisionPreprocessing.sh")
else()
  message(FATAL_ERROR "The current OS is not yet supported")
endif()

add_subdirectory(SRC_DOUBLE)
add_subdirectory(SRC_DOUBLECOMPLEX)
add_subdirectory(SRC_SINGLE)
add_subdirectory(SRC_COMPLEX)
add_subdirectory(EXAMPLE)

# Documentation
if(enable_doc)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/doc/doxygen/doxygen.dox.in
      ${CMAKE_BINARY_DIR}/doxygen.dox @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen.dox
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generating API documentation with doxygen" VERBATIM)
  endif()
endif()

configure_file(${ButterflyPACK_SOURCE_DIR}/make.inc.in ${ButterflyPACK_SOURCE_DIR}/make.inc)

# Add pkg-config support
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/butterflypack.pc.in ${CMAKE_CURRENT_BINARY_DIR}/butterflypack.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/butterflypack.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

set(ConfigPackageLocation lib/cmake/ButterflyPACK)
install(EXPORT MYBPACKTargets FILE butterflypack-targets.cmake
  NAMESPACE ButterflyPACK:: DESTINATION ${ConfigPackageLocation})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(butterflypack-config-version.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY AnyNewerVersion)
install(FILES "${CMAKE_BINARY_DIR}/butterflypack-config-version.cmake"
  DESTINATION ${ConfigPackageLocation})

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/butterflypack-config.cmake.in
  ${CMAKE_BINARY_DIR}/butterflypack-config.cmake
  INSTALL_DESTINATION ${ConfigPackageLocation})

install(FILES "${CMAKE_BINARY_DIR}/butterflypack-config.cmake"
  DESTINATION ${ConfigPackageLocation})
