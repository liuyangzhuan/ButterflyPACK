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
    dBPACK_wrapper.h
    ButterflyPACK_config.fi 	dButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

if (enable_parsec)
  list(APPEND headers dparsec_hmat_ptg_bridge.h)
endif ()

set(common_lapack_sources
dLAPACK_dtrevc3mod.f
dLAPACK_ztrevc3mod.f
dLAPACK_strevc3mod.f
dLAPACK_ctrevc3mod.f
dLAPACK_dgeqp3mod.f
dLAPACK_zgeqp3mod.f
dLAPACK_sgeqp3mod.f
dLAPACK_cgeqp3mod.f
dLAPACK_dgetrfmod.f
dLAPACK_dgetrf2mod.f
dLAPACK_zgetrfmod.f
dLAPACK_zgetrf2mod.f
dLAPACK_sgetrfmod.f
dLAPACK_sgetrf2mod.f
dLAPACK_cgetrfmod.f
dLAPACK_cgetrf2mod.f
)

set(common_scalapack_sources
dSCALAPACK_pdgeqpfmod.f
dSCALAPACK_pzgeqpfmod.f
dSCALAPACK_psgeqpfmod.f
dSCALAPACK_pcgeqpfmod.f
dSCALAPACK_pdgetrfmod.f
dSCALAPACK_pzgetrfmod.f
dSCALAPACK_psgetrfmod.f
dSCALAPACK_pcgetrfmod.f
dSCALAPACK_pdgetf2mod.f
dSCALAPACK_pzgetf2mod.f
dSCALAPACK_psgetf2mod.f
dSCALAPACK_pcgetf2mod.f
)

set(sources
dMAGMA_offset.c
dMAGMA_utilities.f90
dMISC_linkedlist.f90
dMISC_DenseLA.f90
dBPACK_defs.f90
dMISC_utilities.f90
dBplus_utilities.f90
dBPACK_utilities.f90
dBplus_pack_unpack_for_MPI.f90
dBPACK_structure.f90
dBplus_randomized.f90
dBplus_compress.f90
dBplus_factor.f90
dBPACK_solve_mul.f90
dBPACK_constr.f90
dBPACK_randomized.f90
dBPACK_factor.f90
dBPACK_wrapper.f90
dC_BPACK_wrapper.cpp
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
       dparsec_hmat_ptg_bridge.c
       dparsec_hmat_ptg_bridge_fortran.f90)
endif ()


add_library(dbutterflypack OBJECT ${sources} ${headers})

if (enable_parsec)
  set(PARSEC_HMAT_PTG_JDF ${CMAKE_CURRENT_SOURCE_DIR}/dparsec_hmat_factorization_ptg.jdf)
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
  add_custom_target(dbutterflypack_hmat_ptgpp
    DEPENDS ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  add_dependencies(dbutterflypack dbutterflypack_hmat_ptgpp)
  target_sources(dbutterflypack PRIVATE ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  target_compile_definitions(dbutterflypack PRIVATE HAVE_PARSEC)
  target_include_directories(dbutterflypack PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endif ()

if (enable_mpi)
target_link_libraries(dbutterflypack
                       ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
  target_link_libraries(dbutterflypack mpi_dummy ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()

if (enable_parsec)
  target_link_libraries(dbutterflypack PaRSEC::parsec)
endif ()
                      set_target_properties(dbutterflypack PROPERTIES
                      POSITION_INDEPENDENT_CODE ON COMPILE_FLAGS "-DDAT=1"
)

target_include_directories(dbutterflypack PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# ---- Python extension module ----
if (enable_python)
  add_library(dbutterflypack_python SHARED dPy_BPACK_wrapper.cpp)
  set_target_properties(dbutterflypack_python PROPERTIES PREFIX "")
  add_dependencies(dbutterflypack_python dbutterflypack)
  target_link_libraries(dbutterflypack_python
    PUBLIC
      dbutterflypack
      pybind11::module
  )
  target_include_directories(dbutterflypack_python PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
endif()

# Define GNU standard installation directories
include(GNUInstallDirs)

if (enable_python)
  install(TARGETS dbutterflypack_python
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
      dPy_BPACK_wrapper.py
      dPy_BPACK_worker.py
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
  add_custom_target(d_py_bpack ALL
      DEPENDS ${PY_SCRIPT_OUTPUTS}
  )
  install(
      FILES ${PY_SCRIPTS}
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()
