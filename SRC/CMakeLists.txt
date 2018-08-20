set(headers
    C_HODLR_C_Interface.h
    ${PROJECT_BINARY_DIR}/hodlrbf_config.h
)

set(sources
dgeqp3mod.f
zgeqp3mod.f
SCALAPACK_pdgeqpfmod.f
SCALAPACK_pzgeqpfmod.f
Module_file.f90
misc.f90
Utilities.f90
HODLR_structure.f90
Bplus_randomized.f90
Bplus_compress_forward.f90
Bplus_rightmultiply_inverse_randomized.f90
Bplus_inverse_diagonal_randomized_schur_partitioned.f90
HODLR_solve_mul.f90
HODLR_construction.f90
HODLR_randomized.f90
HODLR_factor.f90
HODLR_C_Interface.f90
)



add_library(hodlrbf ${sources} ${HEADERS})
target_link_libraries(hodlrbf
                      ${MPI_Fortran_LIBRARIES} ${LAPACK_LIB} ${SCALAPACK_LIB} m)
set_target_properties(hodlrbf PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR}
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
