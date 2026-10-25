######################################################################
#
# CMakeLists.txt for ButterflyPACK
#
######################################################################

# Required version
cmake_minimum_required(VERSION 3.14...4.3)

# SET(CMAKE_INSTALL_PREFIX < install_path >)

# Project version numbers
project(ButterflyPACK C CXX Fortran)
set(VERSION_MAJOR "5")
set(VERSION_MINOR "0")
set(VERSION_BugFix "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

# Set up options
option(enable_doc       "Build doxygen documentation" OFF)
option(enable_openmp    "Build with OpenMP support" ON)
option(enable_toplevel_openmp    "Use OpenMP at top levels of the solver" OFF)
option(enable_mpi    "Build with MPI support" ON)
option(enable_python "Enable Python" OFF)
option(enable_fftw "Build FFTW" OFF)

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
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
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

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "-DNDEBUG ${CMAKE_Fortran_FLAGS}")
endif()

check_fortran_compiler_flag("-ffree-line-length-none" COMPILER_GNU)
if(COMPILER_GNU)
  set(CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none -ffixed-line-length-none -fimplicit-none -lpthread ${CMAKE_Fortran_FLAGS}")
  check_fortran_compiler_flag("-fallow-argument-mismatch" GNU10)
  if(GNU10)
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch ${CMAKE_Fortran_FLAGS}")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "-frange-check -fbacktrace -fbounds-check -Wconversion ${CMAKE_Fortran_FLAGS}")
  else()
    set(CMAKE_Fortran_FLAGS "-fno-range-check ${CMAKE_Fortran_FLAGS}")
  endif()
endif()

check_fortran_compiler_flag("-no-prec-div" COMPILER_Ifort)
check_fortran_compiler_flag("-qopenmp" COMPILER_Ifx)
if(COMPILER_Ifort OR COMPILER_Ifx)
  set(CMAKE_Fortran_FLAGS "-extend_source 132 -cpp -DIntel ${CMAKE_Fortran_FLAGS}")
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

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
check_cxx_compiler_flag("-qopt-matmul" qoptmatmul)
if(qoptmatmul)
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS  "-qopt-matmul ${CMAKE_CXX_FLAGS}")
  endif()
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

    if(enable_toplevel_openmp)
        set(CMAKE_Fortran_FLAGS "-DHAVE_TOPLEVEL_OPENMP ${OpenMP_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS}")
        set(OpenMP_Fortran_FLAGS_EXPORT "${OpenMP_Fortran_FLAGS}")
        set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ${OpenMP_Fortran_FLAGS}")
    endif()

  endif()




  if(OPENMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "-DHAVE_OPENMP ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(OpenMP_CXX_FLAGS_EXPORT "${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} ${OpenMP_CXX_FLAGS}")
  endif()


# Test OPENMP TASKLOOP (Version 5.0)
  if(COMPILER_PGI OR COMPILER_CRAY OR CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")  # PGI and Cray taskloop seems buggy
  else()
    check_fortran_source_compiles(
"PROGRAM TEST
integer i,sum
sum=0
!$omp parallel
!$omp single
!$omp taskloop default(shared) private(i) reduction(+:sum)
do i=1,10
write(*,*)'iter:',i
sum=sum+i
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
  string(TOLOWER "${BLAS_LIB_STR}" BLAS_LIB_STR_LOWER)
  string(FIND "${BLAS_LIB_STR_LOWER}" "openblas" BUTTERFLYPACK_USE_OPENBLAS)
  if(BUTTERFLYPACK_USE_OPENBLAS GREATER -1)
    set(CMAKE_Fortran_FLAGS "-DHAVE_OPENBLAS_BATCH ${CMAKE_Fortran_FLAGS}")
    message("-- Enabling opt-in OpenBLAS batch GEMM")
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
if(enable_mpi)
  set(CMAKE_Fortran_FLAGS "-DHAVE_MPI ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_CXX_FLAGS "-DHAVE_MPI ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "-DHAVE_MPI ${CMAKE_C_FLAGS}")

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

#--------------------- H2 GPU backend ---------------------
# GPU backend of the format-7 H2 Color factorization (C++ only; independent of
# HAVE_MAGMA, which selects MAGMA kernels in the Fortran butterfly code).
option(enable_h2_gpu "Build the GPU backend of the H2 Color factorization" OFF)
if(enable_h2_gpu)
  if(NOT TPL_H2_MAGMA_LIBRARIES OR NOT TPL_H2_MAGMA_INCLUDE_DIRS)
    message(FATAL_ERROR "enable_h2_gpu requires TPL_H2_MAGMA_LIBRARIES and TPL_H2_MAGMA_INCLUDE_DIRS")
  endif()
  find_package(CUDAToolkit REQUIRED)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 80)  # A100
  endif()
  enable_language(CUDA)
  message("-- Using H2 GPU backend: MAGMA '${TPL_H2_MAGMA_LIBRARIES}', CUDA ${CUDAToolkit_VERSION}, arch ${CMAKE_CUDA_ARCHITECTURES}")
  set(CMAKE_CXX_FLAGS "-DH2_HAVE_GPU -DADD_ ${CMAKE_CXX_FLAGS}")
  include_directories(${TPL_H2_MAGMA_INCLUDE_DIRS} ${CUDAToolkit_INCLUDE_DIRS})
  set(H2_GPU_LIB ${TPL_H2_MAGMA_LIBRARIES} CUDA::cudart CUDA::cusolver CUDA::nvrtc CUDA::cuda_driver)
  #--------------------- nvCOMP (optional) ---------------------
  # Lossy Bitcomp compression of the H2 factors (X_NR is 86-88% of the factor
  # bytes).  Optional: without it the factors are stored uncompressed, as
  # before.  Both libraries are linked, so the host paths
  # (bitcompHost*) are available as well as the device ones.
  if(TPL_H2_NVCOMP_LIBRARIES AND TPL_H2_NVCOMP_INCLUDE_DIRS)
    message("-- Using nvCOMP for H2 factor compression: '${TPL_H2_NVCOMP_LIBRARIES}'")
    set(CMAKE_CXX_FLAGS "-DH2_HAVE_NVCOMP ${CMAKE_CXX_FLAGS}")
    include_directories(${TPL_H2_NVCOMP_INCLUDE_DIRS})
    list(APPEND H2_GPU_LIB ${TPL_H2_NVCOMP_LIBRARIES})
  else()
    message("-- H2 factor compression off (set TPL_H2_NVCOMP_LIBRARIES and TPL_H2_NVCOMP_INCLUDE_DIRS to enable)")
  endif()
endif()

#--------------------- ARPACK ---------------------
if(enable_mpi)
if(TPL_ARPACK_LIBRARIES)
  set(CMAKE_Fortran_FLAGS "-DHAVE_ARPACK ${CMAKE_Fortran_FLAGS}")
  message("-- Using TPL_ARPACK_LIBRARIES='${TPL_ARPACK_LIBRARIES}'")
  set(ARPACK_LIB ${TPL_ARPACK_LIBRARIES})
  # fix up ARPACK library name
  string(REPLACE ";" " " ARPACK_LIB_STR "${ARPACK_LIB}")
  set(ARPACK_LIB_EXPORT ${ARPACK_LIB_STR})
endif()
endif()


#--------------------- ZFP ---------------------
# list(APPEND CMAKE_PREFIX_PATH ${TPL_ZFP_PREFIX}
#     $ENV{ZFP_ROOT} $ENV{ZFP_DIR})
# find_package(ZFP)
# if(ZFP_FOUND)

# set(ZFP_LIB zfp::zfp)
# string (REPLACE ";" " " ZFP_LIB_STR "${ZFP_LIB}")
# set(ZFP_LIB_EXPORT ${ZFP_LIB_STR}) 


# message("ZFP_LIB_STR='${ZFP_LIB_STR}'")

# endif()

if(TPL_ZFP_LIBRARIES)
  set(ZFP_FOUND TRUE)
endif()

if(ZFP_FOUND)
  message("-- Using TPL_ZFP_LIBRARIES='${TPL_ZFP_LIBRARIES}'")
  set(ZFP_LIB ${TPL_ZFP_LIBRARIES})
  # fix up ZFP library name
  string(REPLACE ";" " " ZFP_LIB_STR "${ZFP_LIB}")
  set(ZFP_LIB_EXPORT ${ZFP_LIB_STR})
  set(CMAKE_Fortran_FLAGS "-DHAVE_ZFP ${CMAKE_Fortran_FLAGS} -I${TPL_ZFP_INCLUDE}")
  set(CMAKE_C_FLAGS "-DHAVE_ZFP ${CMAKE_C_FLAGS} -I${TPL_ZFP_INCLUDE}")
endif()

if(enable_fftw)
  if(TPL_FFTW_LIBRARIES)
    message("-- Using TPL_FFTW_LIBRARIES='${TPL_FFTW_LIBRARIES}'")
    set(FFTW_LIB ${TPL_FFTW_LIBRARIES})
  else()
    find_package(FFTW REQUIRED)
    if(FFTW_FOUND)
      set(TPL_FFTW_LIBRARIES "${FFTW_LIBRARIES}" CACHE FILEPATH
      "Set from FindFFTW.cmake FFTW_LIBRARIES." FORCE)
      set(FFTW_LIB ${TPL_FFTW_LIBRARIES})
    endif()
  endif()
  string(REPLACE ";" " " FFTW_LIB_STR "${FFTW_LIB}")
  set(FFTW_LIB_EXPORT ${FFTW_LIB_STR})
  set(CMAKE_Fortran_FLAGS "-DHAVE_FFTW ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "-DHAVE_FFTW ${CMAKE_C_FLAGS}")
endif()


#--------------------- MPI ---------------------
if(enable_mpi)
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
endif()

if(enable_python)
  #--------------------- Python ---------------------
  # Find Python
  find_package(Python REQUIRED COMPONENTS Interpreter Development)
  # Find pybind11
  execute_process(
      COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
      OUTPUT_VARIABLE pybind11_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  list(APPEND CMAKE_PREFIX_PATH "${pybind11_DIR}")
  find_package(pybind11 REQUIRED)
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
if(enable_mpi)
include_directories(${MPI_Fortran_INCLUDE_PATH})
endif()
include_directories(${TPL_LAPACK95_INCLUDE_DIRS})
# bpack_gpu.h, included by the C wrapper headers (BPACK_wrapper.h)
include_directories(${CMAKE_SOURCE_DIR}/GPU_INTERFACE)
# the H2 code (h2_parallel/: core/, structured/, unstructured/) and the GPU
# backends (GPU_BACKEND/: gpu_common/, h2_gpu/, hodlr_gpu/), included by path
# from these roots (e.g. "core/tree.hpp", "gpu_common/device_heap.hpp")
include_directories(${CMAKE_SOURCE_DIR}/h2_parallel ${CMAKE_SOURCE_DIR}/GPU_BACKEND)

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
elseif(UNIX)
  set(USING_MACOS_SED FALSE)
else()
  message(FATAL_ERROR "The current OS is not yet supported")
endif()

# Detect version of SED_TOOL version
message(STATUS "detecting sed version...")
execute_process(
  COMMAND ${SED_TOOL} --version
  WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
  OUTPUT_VARIABLE SED_TOOL_VERSION
  ERROR_QUIET
)
message(STATUS "sed version: ${SED_TOOL_VERSION}")

if(SED_TOOL_VERSION MATCHES ".*GNU sed.*")
    message(STATUS "using GNU sed...")
    set(SED_ARG "") # GNU sed
    set(USING_MACOS_SED FALSE)
endif()

message(STATUS "sed command [${SED_TOOL} -i ${SED_ARG}]")


if(USING_MACOS_SED)
  execute_process(
    COMMAND cp PrecisionPreprocessing.sh PrecisionPreprocessing_mac.sh
    WORKING_DIRECTORY ${ButterflyPACK_SOURCE_DIR}
  )
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


if(NOT enable_mpi)
add_subdirectory(MPI_DUMMY)
endif()
add_subdirectory(SRC_DOUBLE)
add_subdirectory(SRC_DOUBLECOMPLEX)
add_subdirectory(SRC_SINGLE)
add_subdirectory(SRC_COMPLEX)
add_subdirectory(CPP_INTERFACE)

include(GNUInstallDirs)

set(H2_GPU_OBJECTS "")
if (enable_h2_gpu)
  # CUDA kernels of the H2 Color GPU backend (the rest is header-only C++)
  # The text of the entry kernels (GPU_INTERFACE/bpack_gpu_kernels.cuh), which
  # evaluator.cu compiles with NVRTC for entry evaluators given as source
  # text; regenerated when the header changes.
  file(READ ${CMAKE_SOURCE_DIR}/GPU_INTERFACE/bpack_gpu_kernels.cuh BPACK_GPU_KERNELS_TEXT)
  configure_file(${CMAKE_SOURCE_DIR}/GPU_INTERFACE/bpack_gpu_kernels_source.h.in
                 ${CMAKE_BINARY_DIR}/GPU_INTERFACE/bpack_gpu_kernels_source.h @ONLY)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
               ${CMAKE_SOURCE_DIR}/GPU_INTERFACE/bpack_gpu_kernels.cuh)
  add_library(h2_gpu_kernels OBJECT GPU_BACKEND/h2_gpu/device_kernels.cu
                                    GPU_BACKEND/h2_gpu/device_id.cu
                                    GPU_BACKEND/h2_gpu/device_adaptive.cu
                                    GPU_BACKEND/h2_gpu/device_gemm.cu
                                    GPU_BACKEND/h2_gpu/device_sketch.cu
                                    GPU_BACKEND/h2_gpu/solve_kernels.cu
                                    GPU_BACKEND/h2_gpu/h2_matvec_kernels.cu
                                    GPU_BACKEND/h2_gpu/evaluator.cu)
  set_target_properties(h2_gpu_kernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(h2_gpu_kernels PRIVATE H2_HAVE_GPU)
  target_include_directories(h2_gpu_kernels PRIVATE ${CMAKE_BINARY_DIR}/GPU_INTERFACE)
  # --fmad=false (no fused multiply-add contraction) for the device ID only
  # (device_id.cu), needed only by GPU H2 CA0 (replicated CA levels,
  # H2_CA_owner_component 0): a ghost box is eliminated again on each rank
  # that holds it and must give bitwise the same skeleton and factors as its
  # owner's copy, while the device ID picks its launch variant (256 or 1024
  # threads, cooperative) per batch.  Whether a*b+c is contracted to an FMA
  # (one rounding instead of two) is the compiler's choice per template
  # instance, so the variants could round differently; without contraction
  # they give the same bits.  The flag applies to the whole file, so the
  # device ID of Color levels also runs without FMA (a small cost, no
  # accuracy effect); other GPU kernels and CPU code are not affected.
  # The kernels of the adaptive ID rows (device_adaptive.cu, H2_ID_proxy 2)
  # feed the same decisions, so they are compiled the same way.
  set_source_files_properties(GPU_BACKEND/h2_gpu/device_id.cu
                              GPU_BACKEND/h2_gpu/device_adaptive.cu PROPERTIES
                              COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>")
  # CUDA kernels of the HODLR GPU backend (GPU_BACKEND/hodlr_gpu/; the rest is header-only
  # C++ used by SRC/BPACK_gpu_wrapper.cpp)
  add_library(hodlr_gpu_kernels OBJECT GPU_BACKEND/hodlr_gpu/hodlr_kernels.cu)
  set_target_properties(hodlr_gpu_kernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(hodlr_gpu_kernels PRIVATE H2_HAVE_GPU)
  target_include_directories(hodlr_gpu_kernels PRIVATE ${CMAKE_SOURCE_DIR}/GPU_BACKEND)
  # CUDA kernels of the butterfly GPU backend (GPU_BACKEND/bf_gpu/; the rest is
  # header-only C++ used by SRC/BPACK_gpu_wrapper.cpp).  bf_exchange.cu mirrors
  # the sweep's inter-level exchanges and so needs the MPI headers, which the
  # Cray compiler wrapper gives the C++ sources but not nvcc.
  add_library(bf_gpu_kernels OBJECT GPU_BACKEND/bf_gpu/bf_id.cu
                                    GPU_BACKEND/bf_gpu/bf_exchange.cu
                                    GPU_BACKEND/bf_gpu/bf_extract.cu
                                    GPU_BACKEND/bf_gpu/bf_randomized.cu)
  set_target_properties(bf_gpu_kernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(bf_gpu_kernels PRIVATE H2_HAVE_GPU)
  target_include_directories(bf_gpu_kernels PRIVATE ${CMAKE_SOURCE_DIR}/GPU_BACKEND)
  # nvcc does not inherit the compiler wrapper's implicit include path, and
  # these sources reach mpi.h through gpu_runtime.hpp, so it is located once.
  # FindMPI reports no include directory when the wrapper has it built in
  # (Cray), hence the search.
  find_path(BPACK_MPI_INCLUDE_DIR mpi.h
            HINTS ${MPI_CXX_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS} $ENV{MPICH_DIR} $ENV{MPI_ROOT}
                  $ENV{I_MPI_ROOT} $ENV{OPENMPI_DIR}
            PATH_SUFFIXES include)
  if(BPACK_MPI_INCLUDE_DIR)
    message(STATUS "Using mpi.h for the butterfly GPU kernels from '${BPACK_MPI_INCLUDE_DIR}'")
    target_include_directories(bf_gpu_kernels PRIVATE ${BPACK_MPI_INCLUDE_DIR})
  else()
    message(WARNING "mpi.h not found for the butterfly GPU kernels; set BPACK_MPI_INCLUDE_DIR")
  endif()
  set(H2_GPU_OBJECTS $<TARGET_OBJECTS:h2_gpu_kernels> $<TARGET_OBJECTS:hodlr_gpu_kernels>
                     $<TARGET_OBJECTS:bf_gpu_kernels>)
endif()

add_library(butterflypack
  $<TARGET_OBJECTS:dbutterflypack>
  $<TARGET_OBJECTS:zbutterflypack>
  $<TARGET_OBJECTS:sbutterflypack>
  $<TARGET_OBJECTS:cbutterflypack>
  $<TARGET_OBJECTS:butterflypack_cpp>
  ${H2_GPU_OBJECTS})

if (enable_mpi)
  target_link_libraries(butterflypack PUBLIC
    ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES}
    ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${FFTW_LIB} ${SCALAPACK_LIB}
    ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
  target_link_libraries(butterflypack PUBLIC
    mpi_dummy ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${FFTW_LIB}
    ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()

if (enable_parsec)
  find_package(PaRSEC CONFIG REQUIRED)
  target_link_libraries(butterflypack PUBLIC PaRSEC::parsec)
endif()

if (enable_h2_gpu)
  target_link_libraries(butterflypack PUBLIC ${H2_GPU_LIB})
  # Batch-independence test of the GPU kernels (GPU_BACKEND/h2_gpu/
  # GPU_CA_PLAN.md, M1.1); not built by default: make h2_gpu_batch_determinism
  add_executable(h2_gpu_batch_determinism EXCLUDE_FROM_ALL
                 GPU_BACKEND/h2_gpu/tests/batch_determinism.cpp)
  target_compile_definitions(h2_gpu_batch_determinism PRIVATE H2_HAVE_GPU)
  target_link_libraries(h2_gpu_batch_determinism butterflypack)
  # Column-independence test of the device solve's per-box kernels (a column
  # of a many-column launch against that column alone); not built by default:
  # make h2_gpu_solve_columns
  add_executable(h2_gpu_solve_columns EXCLUDE_FROM_ALL
                 GPU_BACKEND/h2_gpu/tests/solve_columns.cpp)
  target_compile_definitions(h2_gpu_solve_columns PRIVATE H2_HAVE_GPU)
  target_link_libraries(h2_gpu_solve_columns butterflypack)
  # ... and of the device H2 matvec's kernels: make h2_gpu_matvec_columns
  add_executable(h2_gpu_matvec_columns EXCLUDE_FROM_ALL
                 GPU_BACKEND/h2_gpu/tests/matvec_columns.cpp)
  target_compile_definitions(h2_gpu_matvec_columns PRIVATE H2_HAVE_GPU)
  target_link_libraries(h2_gpu_matvec_columns butterflypack)
endif()
# Self-test of the CA halo geometry (h2_parallel/core/halo_geometry.hpp, the
# ghost and assisting sets by region against their enumeration; GPU_BACKEND/
# h2_gpu/GPU_CA_PLAN.md, M7b); not built by default: make h2_halo_geometry_test
add_executable(h2_halo_geometry_test EXCLUDE_FROM_ALL h2_parallel/tests/halo_geometry_test.cpp)
target_include_directories(h2_halo_geometry_test PRIVATE ${CMAKE_SOURCE_DIR}/h2_parallel/core)
# Test and timing of the component-owner schedule's graph of the process grid (the local
# construction against the global one; h2_gpu/GPU_CA_PLAN.md, M5e); not built by default: make h2_owner_graph_test
add_executable(h2_owner_graph_test EXCLUDE_FROM_ALL h2_parallel/tests/owner_graph_test.cpp)
target_include_directories(h2_owner_graph_test PRIVATE ${CMAKE_SOURCE_DIR}/h2_parallel/core)
# Self-test of the regions of any number of ranks (h2_parallel/core/region_partition.hpp, the
# unstructured Color path; h2_gpu/GPU_CA_PLAN.md); not built by default: make h2_region_partition_test
add_executable(h2_region_partition_test EXCLUDE_FROM_ALL h2_parallel/tests/region_partition_test.cpp)
target_include_directories(h2_region_partition_test PRIVATE ${CMAKE_SOURCE_DIR}/h2_parallel/core)

set_target_properties(butterflypack PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${VERSION_MAJOR})

target_include_directories(butterflypack PUBLIC
  $<BUILD_INTERFACE:${ButterflyPACK_SOURCE_DIR}/GPU_INTERFACE>
  $<BUILD_INTERFACE:${ButterflyPACK_SOURCE_DIR}/CPP_INTERFACE>
  $<BUILD_INTERFACE:${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLE>
  $<BUILD_INTERFACE:${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLECOMPLEX>
  $<BUILD_INTERFACE:${ButterflyPACK_SOURCE_DIR}/SRC_SINGLE>
  $<BUILD_INTERFACE:${ButterflyPACK_SOURCE_DIR}/SRC_COMPLEX>
  $<BUILD_INTERFACE:${ButterflyPACK_BINARY_DIR}/CPP_INTERFACE>
  $<BUILD_INTERFACE:${ButterflyPACK_BINARY_DIR}/SRC_DOUBLE>
  $<BUILD_INTERFACE:${ButterflyPACK_BINARY_DIR}/SRC_DOUBLECOMPLEX>
  $<BUILD_INTERFACE:${ButterflyPACK_BINARY_DIR}/SRC_SINGLE>
  $<BUILD_INTERFACE:${ButterflyPACK_BINARY_DIR}/SRC_COMPLEX>
  $<INSTALL_INTERFACE:include>)

install(TARGETS butterflypack EXPORT MYBPACKTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# the interface of the GPU evaluators (doc/gpu_kernels.md)
install(FILES GPU_INTERFACE/bpack_gpu.h GPU_INTERFACE/bpack_gpu_entry.cuh GPU_INTERFACE/bpack_gpu_kernels.cuh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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

set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/ButterflyPACK)
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
