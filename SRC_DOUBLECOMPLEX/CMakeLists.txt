set(headers
    zC_HODLR_wrapper.h
    # HODLR_config.fi
    # ${PROJECT_BINARY_DIR}/HODLR_config.h
)

set(sources
../EXAMPLE/EMCURV_Module.f90
../EXAMPLE/EMSURF_Module.f90
zdgeqp3mod.f
zzgeqp3mod.f
zSCALAPACK_pdgeqpfmod.f
zSCALAPACK_pzgeqpfmod.f
zDenseLA.f90
zHODLR_defs.f90
zmisc.f90
zHODLR_utilities.f90
zHODLR_structure.f90
zBplus_randomized.f90
zBplus_compress.f90
zBplus_update.f90
zBplus_inversion.f90
zHODLR_solve_mul.f90
zHODLR_constr.f90
zHODLR_randomized.f90
zHODLR_factor.f90
zHODLR_wrapper.f90
)



add_library(zhodlrbf ${sources} ${HEADERS})
target_link_libraries(zhodlrbf
                      ${MPI_Fortran_LIBRARIES} ${LAPACK_LIB} ${SCALAPACK_LIB} m)
set_target_properties(zhodlrbf PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} COMPILE_FLAGS "-DDAT=0"
)


# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS zhodlrbf
# DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
     DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${headers}
# DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
