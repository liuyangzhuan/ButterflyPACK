set(headers
    C_BPACK_wrapper.h
    ButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

set(sources
MAGMA_offset.c
MAGMA_utilities.f90
LAPACK_dgeqp3mod.f
LAPACK_zgeqp3mod.f
LAPACK_sgeqp3mod.f
LAPACK_cgeqp3mod.f
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
)

if (enable_mpi)
  list(APPEND sources SCALAPACK_pdgeqpfmod.f SCALAPACK_pzgeqpfmod.f SCALAPACK_psgeqpfmod.f SCALAPACK_pcgeqpfmod.f)
endif ()


add_library(butterflypack ${sources} ${HEADERS})
target_link_libraries(butterflypack
                      ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
set_target_properties(butterflypack PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} COMPILE_FLAGS "-DDAT"
)

target_include_directories(butterflypack PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS butterflypack EXPORT MYBPACKTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


#SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.mod")