include_directories(${ButterflyPACK_SOURCE_DIR}/MPI_DUMMY)

option(enable_parsec "Build optional PaRSEC PTG support" OFF)

if (enable_parsec)
  if (NOT enable_mpi)
    message(FATAL_ERROR "enable_parsec requires enable_mpi")
  endif ()
  find_package(PaRSEC CONFIG REQUIRED)
  if (NOT TARGET PaRSEC::parsec)
    message(FATAL_ERROR "enable_parsec requires the PaRSEC::parsec CMake target")
  endif ()
  if (NOT TARGET PaRSEC::parsec-ptgpp)
    message(FATAL_ERROR "enable_parsec requires the PaRSEC::parsec-ptgpp CMake target")
  endif ()
endif ()

set(headers
    sBPACK_wrapper.h
    ButterflyPACK_config.fi 	sButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

if (enable_parsec)
  list(APPEND headers sparsec_hmat_ptg_bridge.h)
endif ()

set(common_lapack_sources
sLAPACK_dtrevc3mod.f
sLAPACK_ztrevc3mod.f
sLAPACK_strevc3mod.f
sLAPACK_ctrevc3mod.f
sLAPACK_dgeqp3mod.f
sLAPACK_zgeqp3mod.f
sLAPACK_sgeqp3mod.f
sLAPACK_cgeqp3mod.f
sLAPACK_dgetrfmod.f
sLAPACK_dgetrf2mod.f
sLAPACK_zgetrfmod.f
sLAPACK_zgetrf2mod.f
sLAPACK_sgetrfmod.f
sLAPACK_sgetrf2mod.f
sLAPACK_cgetrfmod.f
sLAPACK_cgetrf2mod.f
)

set(common_scalapack_sources
sSCALAPACK_pdgeqpfmod.f
sSCALAPACK_pzgeqpfmod.f
sSCALAPACK_psgeqpfmod.f
sSCALAPACK_pcgeqpfmod.f
sSCALAPACK_pdgetrfmod.f
sSCALAPACK_pzgetrfmod.f
sSCALAPACK_psgetrfmod.f
sSCALAPACK_pcgetrfmod.f
sSCALAPACK_pdgetf2mod.f
sSCALAPACK_pzgetf2mod.f
sSCALAPACK_psgetf2mod.f
sSCALAPACK_pcgetf2mod.f
)

set(sources
sMAGMA_offset.c
sMAGMA_utilities.f90
sMISC_linkedlist.f90
sMISC_DenseLA.f90
sBPACK_defs.f90
sMISC_utilities.f90
sBplus_utilities.f90
sBPACK_utilities.f90
sBplus_pack_unpack_for_MPI.f90
sBPACK_structure.f90
sBplus_randomized.f90
sBplus_compress.f90
sBplus_factor.f90
sBPACK_solve_mul.f90
sBPACK_constr.f90
sBPACK_randomized.f90
sBPACK_factor.f90
sBPACK_wrapper.f90
sC_BPACK_wrapper.cpp
)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "/SRC$" OR
    CMAKE_CURRENT_SOURCE_DIR MATCHES "/SRC_DOUBLE$")
  list(APPEND sources ${common_lapack_sources})
  if (enable_mpi)
    list(APPEND sources ${common_scalapack_sources})
  endif ()
endif ()

if (enable_parsec)
  list(APPEND sources
       sparsec_hmat_ptg_bridge.c
       sparsec_hmat_ptg_bridge_fortran.f90)
endif ()


add_library(sbutterflypack OBJECT ${sources} ${headers})

if (enable_parsec)
  set(PARSEC_HMAT_PTG_JDF ${CMAKE_CURRENT_SOURCE_DIR}/sparsec_hmat_factorization_ptg.jdf)
  get_filename_component(PARSEC_HMAT_PTG_NAME ${PARSEC_HMAT_PTG_JDF} NAME_WE)
  set(PARSEC_HMAT_PTG_C ${CMAKE_CURRENT_BINARY_DIR}/${PARSEC_HMAT_PTG_NAME}.c)
  set(PARSEC_HMAT_PTG_H ${CMAKE_CURRENT_BINARY_DIR}/${PARSEC_HMAT_PTG_NAME}.h)

  add_custom_command(
    OUTPUT ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C}
    COMMAND $<TARGET_FILE:PaRSEC::parsec-ptgpp>
            --noline --Wremoteref -E
            -i ${PARSEC_HMAT_PTG_JDF}
            -C ${PARSEC_HMAT_PTG_C}
            -H ${PARSEC_HMAT_PTG_H}
            -f ${PARSEC_HMAT_PTG_NAME}
    MAIN_DEPENDENCY ${PARSEC_HMAT_PTG_JDF}
    DEPENDS ${PARSEC_HMAT_PTG_JDF} PaRSEC::parsec-ptgpp
    VERBATIM)
  add_custom_target(sbutterflypack_hmat_ptgpp
    DEPENDS ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  add_dependencies(sbutterflypack sbutterflypack_hmat_ptgpp)
  target_sources(sbutterflypack PRIVATE ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  target_compile_definitions(sbutterflypack PRIVATE HAVE_PARSEC)
  target_include_directories(sbutterflypack PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endif ()

if (enable_mpi)
target_link_libraries(sbutterflypack
                       ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
  target_link_libraries(sbutterflypack mpi_dummy ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()

if (enable_parsec)
  target_link_libraries(sbutterflypack PaRSEC::parsec)
endif ()
                      set_target_properties(sbutterflypack PROPERTIES
                      POSITION_INDEPENDENT_CODE ON COMPILE_FLAGS "-DDAT=3"
)

target_include_directories(sbutterflypack PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# ---- Python extension module ----
if (enable_python)
  add_library(sbutterflypack_python SHARED sPy_BPACK_wrapper.cpp)
  set_target_properties(sbutterflypack_python PROPERTIES PREFIX "")
  add_dependencies(sbutterflypack_python sbutterflypack)
  target_link_libraries(sbutterflypack_python
    PUBLIC
      sbutterflypack
      pybind11::module
  )
  target_include_directories(sbutterflypack_python PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
endif()

# Define GNU standard installation directories
include(GNUInstallDirs)

if (enable_python)
  install(TARGETS sbutterflypack_python
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


#SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.mod")


if (enable_python)
  set(PY_SCRIPTS
      sPy_BPACK_wrapper.py
      sPy_BPACK_worker.py
  )
  set(PY_SCRIPT_OUTPUTS)
  foreach(pyfile IN LISTS PY_SCRIPTS)
      set(src ${CMAKE_CURRENT_SOURCE_DIR}/${pyfile})
      set(dst ${CMAKE_CURRENT_BINARY_DIR}/${pyfile})

      add_custom_command(
          OUTPUT ${dst}
          COMMAND ${CMAKE_COMMAND} -E copy_if_different ${src} ${dst}
          DEPENDS ${src}
          COMMENT "Copying ${pyfile}"
      )

      list(APPEND PY_SCRIPT_OUTPUTS ${dst})
  endforeach()
  add_custom_target(s_py_bpack ALL
      DEPENDS ${PY_SCRIPT_OUTPUTS}
  )
  install(
      FILES ${PY_SCRIPTS}
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()
