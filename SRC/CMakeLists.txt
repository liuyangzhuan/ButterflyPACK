set(headers
    C_HODLR_wrapper.h
    # HODLR_config.fi
    # ${PROJECT_BINARY_DIR}/HODLR_config.h
)

set(sources
dgeqp3mod.f
zgeqp3mod.f
SCALAPACK_pdgeqpfmod.f
SCALAPACK_pzgeqpfmod.f
DenseLA.f90
HODLR_defs.f90
misc.f90
HODLR_utilities.f90
HODLR_structure.f90
Bplus_randomized.f90
Bplus_compress.f90
Bplus_update.f90
Bplus_inversion.f90
HODLR_solve_mul.f90
HODLR_constr.f90
HODLR_randomized.f90
HODLR_factor.f90
HODLR_wrapper.f90
)



add_library(hodlrbf ${sources} ${HEADERS})
target_link_libraries(hodlrbf
                      ${MPI_Fortran_LIBRARIES} ${LAPACK_LIB} ${SCALAPACK_LIB} m)
set_target_properties(hodlrbf PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} COMPILE_FLAGS "-DDAT"
)


# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS hodlrbf
# DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
     DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${headers}
# DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
