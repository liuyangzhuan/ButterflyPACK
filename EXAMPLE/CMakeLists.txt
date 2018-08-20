include_directories(${HODLRBF_SOURCE_DIR}/SRC)
include_directories( ${HODLRBF_BINARY_DIR}/SRC)
# Libs linked to all of the examples
set(all_link_libs hodlrbf ${MPI_Fortran_LIBRARIES} ${LAPACK_LIB} ${SCALAPACK_LIB} m)

if (enable_complex)
set(FZEXM1 EMCURV_Driver.f90)
add_executable(fem2d ${FZEXM1})
target_link_libraries(fem2d ${all_link_libs})

set(FZEXM2 EMSURF_Driver.f90)
add_executable(fem3d ${FZEXM2})
target_link_libraries(fem3d ${all_link_libs})

# set(FZEXM3 HODLR_Randconst.f90)
# add_executable(frandh ${FZEXM3})
# target_link_libraries(frandh ${all_link_libs})

else()
set(CDEXM InterfaceTest.cpp)
add_executable(ctest ${CDEXM})
target_link_libraries(ctest ${all_link_libs})

set(FDEXM1 FULLMAT_Driver.f90)
add_executable(ffull ${FDEXM1})
target_link_libraries(ffull ${all_link_libs})

set(FDEXM2 KERREG_Driver.f90)
add_executable(fkerreg ${FDEXM2})
target_link_libraries(fkerreg ${all_link_libs})

endif()





