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
    cBPACK_wrapper.h
    ButterflyPACK_config.fi 	cButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

if (enable_parsec)
  list(APPEND headers cparsec_hmat_ptg_bridge.h)
endif ()

set(common_lapack_sources
cLAPACK_dtrevc3mod.f
cLAPACK_ztrevc3mod.f
cLAPACK_strevc3mod.f
cLAPACK_ctrevc3mod.f
cLAPACK_dgeqp3mod.f
cLAPACK_zgeqp3mod.f
cLAPACK_sgeqp3mod.f
cLAPACK_cgeqp3mod.f
cLAPACK_dgetrfmod.f
cLAPACK_dgetrf2mod.f
cLAPACK_zgetrfmod.f
cLAPACK_zgetrf2mod.f
cLAPACK_sgetrfmod.f
cLAPACK_sgetrf2mod.f
cLAPACK_cgetrfmod.f
cLAPACK_cgetrf2mod.f
)

set(common_scalapack_sources
cSCALAPACK_pdgeqpfmod.f
cSCALAPACK_pzgeqpfmod.f
cSCALAPACK_psgeqpfmod.f
cSCALAPACK_pcgeqpfmod.f
cSCALAPACK_pdgetrfmod.f
cSCALAPACK_pzgetrfmod.f
cSCALAPACK_psgetrfmod.f
cSCALAPACK_pcgetrfmod.f
cSCALAPACK_pdgetf2mod.f
cSCALAPACK_pzgetf2mod.f
cSCALAPACK_psgetf2mod.f
cSCALAPACK_pcgetf2mod.f
)

set(sources
cMAGMA_offset.c
cMAGMA_utilities.f90
cMISC_linkedlist.f90
cMISC_DenseLA.f90
cBPACK_defs.f90
cMISC_utilities.f90
cBplus_utilities.f90
cBPACK_utilities.f90
cBplus_pack_unpack_for_MPI.f90
cBPACK_structure.f90
cBplus_randomized.f90
cBplus_compress.f90
cBplus_factor.f90
cBPACK_solve_mul.f90
cBPACK_constr.f90
cBPACK_randomized.f90
cBPACK_factor.f90
cBPACK_wrapper.f90
cC_BPACK_wrapper.cpp
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
       cparsec_hmat_ptg_bridge.c
       cparsec_hmat_ptg_bridge_fortran.f90)
endif ()


add_library(cbutterflypack OBJECT ${sources} ${headers})

if (enable_parsec)
  set(PARSEC_HMAT_PTG_JDF ${CMAKE_CURRENT_SOURCE_DIR}/cparsec_hmat_factorization_ptg.jdf)
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
  add_custom_target(cbutterflypack_hmat_ptgpp
    DEPENDS ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  add_dependencies(cbutterflypack cbutterflypack_hmat_ptgpp)
  target_sources(cbutterflypack PRIVATE ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  target_compile_definitions(cbutterflypack PRIVATE HAVE_PARSEC)
  target_include_directories(cbutterflypack PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endif ()

if (enable_mpi)
target_link_libraries(cbutterflypack
                       ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
  target_link_libraries(cbutterflypack mpi_dummy ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()

if (enable_parsec)
  target_link_libraries(cbutterflypack PaRSEC::parsec)
endif ()
                      set_target_properties(cbutterflypack PROPERTIES
                      POSITION_INDEPENDENT_CODE ON COMPILE_FLAGS "-DDAT=2"
)

target_include_directories(cbutterflypack PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# ---- Python extension module ----
if (enable_python)
  add_library(cbutterflypack_python SHARED cPy_BPACK_wrapper.cpp)
  set_target_properties(cbutterflypack_python PROPERTIES PREFIX "")
  add_dependencies(cbutterflypack_python cbutterflypack)
  target_link_libraries(cbutterflypack_python
    PUBLIC
      cbutterflypack
      pybind11::module
  )
  target_include_directories(cbutterflypack_python PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
endif()

# Define GNU standard installation directories
include(GNUInstallDirs)

if (enable_python)
  install(TARGETS cbutterflypack_python
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
      cPy_BPACK_wrapper.py
      cPy_BPACK_worker.py
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
  add_custom_target(c_py_bpack ALL
      DEPENDS ${PY_SCRIPT_OUTPUTS}
  )
  install(
      FILES ${PY_SCRIPTS}
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()
