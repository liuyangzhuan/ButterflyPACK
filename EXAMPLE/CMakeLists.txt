include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLECOMPLEX)
include_directories(${ButterflyPACK_SOURCE_DIR}/EXAMPLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLECOMPLEX)
# Libs linked to all of the examples
set(all_link_libs_z zbutterflypack ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} m)
set(all_link_libs_d dbutterflypack ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} m)

# if (enable_complex)
set(FZEXM1 EMCURV_Driver.f90 EMCURV_Module.f90)
add_executable(ie2d ${FZEXM1})
target_link_libraries(ie2d ${all_link_libs_z})

set(FZEXM2 EMSURF_Driver.f90 EMSURF_Module.f90)
add_executable(ie3d ${FZEXM2})
target_link_libraries(ie3d ${all_link_libs_z})

set(FZEXM3 SMAT_Driver.f90)
add_executable(smat ${FZEXM3})
target_link_libraries(smat ${all_link_libs_z})

set(FZEXM4 EMCURV_Eigen_Driver.f90 EMCURV_Module.f90)
add_executable(ie2deigen ${FZEXM4})
target_link_libraries(ie2deigen ${all_link_libs_z})

set(FZEXM5 FrontalDist_Driver.f90)
add_executable(frontaldist ${FZEXM5})
target_link_libraries(frontaldist ${all_link_libs_z})

set(FZEXM6 EMSURF_Eigen_Driver.f90 EMSURF_Module.f90)
add_executable(ie3deigen ${FZEXM6})
target_link_libraries(ie3deigen ${all_link_libs_z})

# else()
set(CDEXM InterfaceTest.cpp)
add_executable(ctest ${CDEXM})
target_link_libraries(ctest ${all_link_libs_d})


set(FDEXM1 FULLMAT_Driver.f90)
add_executable(full ${FDEXM1})
target_link_libraries(full ${all_link_libs_d})
target_link_libraries(full ${all_link_libs_z})

set(FDEXM2 KERREG_Driver.f90)
add_executable(krr ${FDEXM2})
target_link_libraries(krr ${all_link_libs_d})

# set(FDEXM3 FULLMATKERREG_Driver.f90)
# add_executable(fullkrr ${FDEXM3})
# target_link_libraries(fullkrr ${all_link_libs_d})

set(FDEXM4 Frontal_Driver.f90)
add_executable(frontal ${FDEXM4})
target_link_libraries(frontal ${all_link_libs_d})


# endif()





