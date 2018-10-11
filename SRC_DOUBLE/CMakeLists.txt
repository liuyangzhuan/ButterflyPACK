set(headers
    dC_HODLR_wrapper.h
    # HODLR_config.fi
    # ${PROJECT_BINARY_DIR}/HODLR_config.h
)

set(sources
ddgeqp3mod.f
dzgeqp3mod.f
dSCALAPACK_pdgeqpfmod.f
dSCALAPACK_pzgeqpfmod.f
dDenseLA.f90
dHODLR_defs.f90
dmisc.f90
dHODLR_utilities.f90
dHODLR_structure.f90
dBplus_randomized.f90
dBplus_compress.f90
dBplus_update.f90
dBplus_inversion.f90
dHODLR_solve_mul.f90
dHODLR_constr.f90
dHODLR_randomized.f90
dHODLR_factor.f90
dHODLR_wrapper.f90
)



add_library(dhodlrbf ${sources} ${HEADERS})
target_link_libraries(dhodlrbf
                      ${MPI_Fortran_LIBRARIES} ${LAPACK_LIB} ${SCALAPACK_LIB} m)
set_target_properties(dhodlrbf PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} COMPILE_FLAGS "-DDAT=1"
)


# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS dhodlrbf
# DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
     DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${headers}
# DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
