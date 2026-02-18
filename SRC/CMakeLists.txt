include_directories(${ButterflyPACK_SOURCE_DIR}/MPI_DUMMY)

set(headers
    BPACK_wrapper.h
    ButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

set(sources
MAGMA_offset.c
MAGMA_utilities.f90
LAPACK_dtrevc3mod.f
LAPACK_ztrevc3mod.f
LAPACK_strevc3mod.f
LAPACK_ctrevc3mod.f
LAPACK_dgeqp3mod.f
LAPACK_zgeqp3mod.f
LAPACK_sgeqp3mod.f
LAPACK_cgeqp3mod.f
LAPACK_dgetrfmod.f
LAPACK_dgetrf2mod.f
LAPACK_zgetrfmod.f
LAPACK_zgetrf2mod.f
LAPACK_sgetrfmod.f
LAPACK_sgetrf2mod.f
LAPACK_cgetrfmod.f
LAPACK_cgetrf2mod.f
MISC_linkedlist.f90
MISC_DenseLA.f90
BPACK_defs.f90
MISC_utilities.f90
Bplus_utilities.f90
BPACK_utilities.f90
Bplus_pack_unpack_for_MPI.f90
BPACK_structure.f90
Bplus_randomized.f90
Bplus_compress.f90
Bplus_factor.f90
BPACK_solve_mul.f90
BPACK_constr.f90
BPACK_randomized.f90
BPACK_factor.f90
BPACK_wrapper.f90
C_BPACK_wrapper.cpp
)

if (enable_mpi)
  list(APPEND sources SCALAPACK_pdgeqpfmod.f SCALAPACK_pzgeqpfmod.f SCALAPACK_psgeqpfmod.f SCALAPACK_pcgeqpfmod.f SCALAPACK_pdgetrfmod.f SCALAPACK_pzgetrfmod.f SCALAPACK_psgetrfmod.f SCALAPACK_pcgetrfmod.f SCALAPACK_pdgetf2mod.f SCALAPACK_pzgetf2mod.f SCALAPACK_psgetf2mod.f SCALAPACK_pcgetf2mod.f)
endif ()


add_library(butterflypack ${sources} ${HEADERS})

if (enable_mpi)
target_link_libraries(butterflypack
                       ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
target_link_libraries(butterflypack mpi_dummy ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()
                      set_target_properties(butterflypack PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} COMPILE_FLAGS "-DDAT"
)

target_include_directories(butterflypack PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)


# ---- Python extension module ----
if (enable_python)
  add_library(butterflypack_python SHARED Py_BPACK_wrapper.cpp)
  set_target_properties(butterflypack_python PROPERTIES PREFIX "")
  add_dependencies(butterflypack_python butterflypack)
  target_link_libraries(butterflypack_python
    PUBLIC
      butterflypack
      pybind11::module
  )
  target_include_directories(butterflypack_python PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
endif()

# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS butterflypack EXPORT MYBPACKTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if (enable_python)
  install(TARGETS butterflypack_python EXPORT MYBPACKTargets
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
      Py_BPACK_wrapper.py
      Py_BPACK_worker.py
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
  add_custom_target(py_bpack ALL
      DEPENDS ${PY_SCRIPT_OUTPUTS}
  )
  install(
      FILES ${PY_SCRIPTS}
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()