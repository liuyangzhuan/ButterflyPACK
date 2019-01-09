set(headers
    C_BPACK_wrapper.h
    # ButterflyPACK_config.fi
    # ${PROJECT_BINARY_DIR}/ButterflyPACK_config.h
)

set(sources
dgeqp3mod.f
zgeqp3mod.f
SCALAPACK_pdgeqpfmod.f
SCALAPACK_pzgeqpfmod.f
DenseLA.f90
BPACK_defs.f90
misc.f90
Bplus_utilities.f90
BPACK_utilities.f90
blocks_pack_unpack_for_MPI.f90
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



add_library(butterflypack ${sources} ${HEADERS})
target_link_libraries(butterflypack
                      ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${SCALAPACK_LIB} m)
set_target_properties(butterflypack PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} COMPILE_FLAGS "-DDAT"
)


# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS butterflypack
# DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
     DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${headers}
# DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
