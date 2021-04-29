include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLECOMPLEX)
include_directories(${ButterflyPACK_SOURCE_DIR}/EXAMPLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLECOMPLEX)
# Libs linked to all of the examples
set(all_link_libs_z zbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB}  m)
set(all_link_libs_d dbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)


########################################### double complex Fortran drivers

set(FZEXM1 EMCURV_Driver.f90 EMCURV_Module.f90)
add_executable(ie2d ${FZEXM1})
target_link_libraries(ie2d ${all_link_libs_z})

set(FZEXM2 EMCURV_Eigen_Driver.f90 EMCURV_Module.f90)
add_executable(ie2deigen ${FZEXM2})
target_link_libraries(ie2deigen ${all_link_libs_z})
add_dependencies(ie2deigen ie2d)


set(FZEXM3 EMSURF_Driver.f90 EMSURF_Module.f90)
add_executable(ie3d ${FZEXM3})
target_link_libraries(ie3d ${all_link_libs_z})

set(FZEXM4 EMSURF_Eigen_Driver.f90 EMSURF_Module.f90)
add_executable(ie3deigen ${FZEXM4})
target_link_libraries(ie3deigen ${all_link_libs_z})
add_dependencies(ie3deigen ie3d)


set(FZEXM5 SMAT_Driver.f90)
add_executable(smat ${FZEXM5})
target_link_libraries(smat ${all_link_libs_z})

set(FZEXM6 FrontalDist_Driver.f90)
add_executable(frontaldist ${FZEXM6})
target_link_libraries(frontaldist ${all_link_libs_z})
add_dependencies(frontaldist smat)


set(FZEXM7 Taylor2D.cpp G2D/G2D.cpp G2D/rtGetInf.cpp G2D/rtGetNaN.cpp G2D/rt_nonfinite.cpp G2D/bessel.cpp)
add_executable(go2d ${FZEXM7})
target_link_libraries(go2d ${all_link_libs_z})
add_dependencies(go2d frontaldist)

set(FZEXM8 Taylor3D.cpp G3D/G3D.cpp G3D/rtGetInf.cpp G3D/rtGetNaN.cpp G3D/rt_nonfinite.cpp)
add_executable(go3d ${FZEXM8})
target_link_libraries(go3d ${all_link_libs_z})
add_dependencies(go3d go2d)


set(FZEXM9 EMSURF_Port_Eigen_Driver.f90 EMSURF_Port_Module.f90)
add_executable(ie3dporteigen ${FZEXM9})
target_link_libraries(ie3dporteigen ${all_link_libs_z})
add_dependencies(ie3dporteigen ie3d)


########################################### double Fortran drivers

set(FDEXM1 FULLMAT_Driver.f90)
add_executable(full ${FDEXM1})
target_link_libraries(full ${all_link_libs_d})
target_link_libraries(full ${all_link_libs_z})
add_dependencies(full frontaldist)

set(FDEXM1-1 FULLMAT_Driver_simple.f90)
add_executable(full_simple ${FDEXM1-1})
target_link_libraries(full_simple ${all_link_libs_d})
target_link_libraries(full_simple ${all_link_libs_z})
add_dependencies(full_simple full)

set(FDEXM2 KERREG_Driver.f90)
add_executable(krr ${FDEXM2})
target_link_libraries(krr ${all_link_libs_d})
add_dependencies(krr full_simple)

set(FDEXM3 Frontal_Driver.f90)
add_executable(frontal ${FDEXM3})
target_link_libraries(frontal ${all_link_libs_d})
add_dependencies(frontal krr)

set(FDEXM4 FULLMATKERREG_Driver.f90)
add_executable(fullkrr ${FDEXM4})
target_link_libraries(fullkrr ${all_link_libs_d})
add_dependencies(fullkrr frontal)



########################################### double C++ drivers

set(CDEXM InterfaceTest.cpp)
add_executable(ctest ${CDEXM})
target_link_libraries(ctest ${all_link_libs_d})

set(CDEXM InterfaceTest_simple.cpp)
add_executable(ctest_simple ${CDEXM})
target_link_libraries(ctest_simple ${all_link_libs_d})
add_dependencies(ctest_simple ctest)

########################################### double complex C++ drivers

set(CFEXM FIO_Driver.cpp)
add_executable(cfio ${CFEXM})
target_link_libraries(cfio ${all_link_libs_d})
target_link_libraries(cfio ${all_link_libs_z})
add_dependencies(cfio ctest)







