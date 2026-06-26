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
    zBPACK_wrapper.h
    ButterflyPACK_config.fi 	zButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

if (enable_parsec)
  list(APPEND headers zparsec_hmat_ptg_bridge.h)
endif ()

set(common_lapack_sources
zLAPACK_dtrevc3mod.f
zLAPACK_ztrevc3mod.f
zLAPACK_strevc3mod.f
zLAPACK_ctrevc3mod.f
zLAPACK_dgeqp3mod.f
zLAPACK_zgeqp3mod.f
zLAPACK_sgeqp3mod.f
zLAPACK_cgeqp3mod.f
zLAPACK_dgetrfmod.f
zLAPACK_dgetrf2mod.f
zLAPACK_zgetrfmod.f
zLAPACK_zgetrf2mod.f
zLAPACK_sgetrfmod.f
zLAPACK_sgetrf2mod.f
zLAPACK_cgetrfmod.f
zLAPACK_cgetrf2mod.f
)

set(common_scalapack_sources
zSCALAPACK_pdgeqpfmod.f
zSCALAPACK_pzgeqpfmod.f
zSCALAPACK_psgeqpfmod.f
zSCALAPACK_pcgeqpfmod.f
zSCALAPACK_pdgetrfmod.f
zSCALAPACK_pzgetrfmod.f
zSCALAPACK_psgetrfmod.f
zSCALAPACK_pcgetrfmod.f
zSCALAPACK_pdgetf2mod.f
zSCALAPACK_pzgetf2mod.f
zSCALAPACK_psgetf2mod.f
zSCALAPACK_pcgetf2mod.f
)

set(sources
zMAGMA_offset.c
zMAGMA_utilities.f90
zMISC_linkedlist.f90
zMISC_DenseLA.f90
zBPACK_defs.f90
zMISC_utilities.f90
zBplus_utilities.f90
zBPACK_utilities.f90
zBplus_pack_unpack_for_MPI.f90
zBPACK_structure.f90
zBplus_randomized.f90
zBplus_compress.f90
zBplus_factor.f90
zBPACK_solve_mul.f90
zBPACK_constr.f90
zBPACK_randomized.f90
zBPACK_factor.f90
zBPACK_wrapper.f90
zC_BPACK_wrapper.cpp
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
       zparsec_hmat_ptg_bridge.c
       zparsec_hmat_ptg_bridge_fortran.f90)
endif ()


add_library(zbutterflypack OBJECT ${sources} ${headers})

if (enable_parsec)
  set(PARSEC_HMAT_PTG_JDF ${CMAKE_CURRENT_SOURCE_DIR}/zparsec_hmat_factorization_ptg.jdf)
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
  add_custom_target(zbutterflypack_hmat_ptgpp
    DEPENDS ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  add_dependencies(zbutterflypack zbutterflypack_hmat_ptgpp)
  target_sources(zbutterflypack PRIVATE ${PARSEC_HMAT_PTG_H} ${PARSEC_HMAT_PTG_C})
  target_compile_definitions(zbutterflypack PRIVATE HAVE_PARSEC)
  target_include_directories(zbutterflypack PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endif ()

if (enable_mpi)
target_link_libraries(zbutterflypack
                       ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
  target_link_libraries(zbutterflypack mpi_dummy ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()

if (enable_parsec)
  target_link_libraries(zbutterflypack PaRSEC::parsec)
endif ()
                      set_target_properties(zbutterflypack PROPERTIES
                      POSITION_INDEPENDENT_CODE ON COMPILE_FLAGS "-DDAT=0"
)

target_include_directories(zbutterflypack PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# ---- Python extension module ----
if (enable_python)
  add_library(zbutterflypack_python SHARED zPy_BPACK_wrapper.cpp)
  set_target_properties(zbutterflypack_python PROPERTIES PREFIX "")
  add_dependencies(zbutterflypack_python zbutterflypack)
  target_link_libraries(zbutterflypack_python
    PUBLIC
      zbutterflypack
      pybind11::module
  )
  target_include_directories(zbutterflypack_python PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
endif()

# Define GNU standard installation directories
include(GNUInstallDirs)

if (enable_python)
  install(TARGETS zbutterflypack_python
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
      zPy_BPACK_wrapper.py
      zPy_BPACK_worker.py
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
  add_custom_target(z_py_bpack ALL
      DEPENDS ${PY_SCRIPT_OUTPUTS}
  )
  install(
      FILES ${PY_SCRIPTS}
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()
