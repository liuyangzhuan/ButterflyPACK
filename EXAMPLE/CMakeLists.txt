include_directories(${ButterflyPACK_SOURCE_DIR}/CPP_INTERFACE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLECOMPLEX)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_SINGLE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_COMPLEX)
include_directories(${ButterflyPACK_SOURCE_DIR}/EXAMPLE)
if(NOT enable_mpi)
include_directories(${ButterflyPACK_SOURCE_DIR}/MPI_DUMMY)
endif()
include_directories( ${ButterflyPACK_BINARY_DIR}/CPP_INTERFACE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLECOMPLEX)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_SINGLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_COMPLEX)
# Libs linked to all of the examples
if(enable_mpi)
set(all_link_libs_z zbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB}  m)
set(all_link_libs_d dbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
set(all_link_libs_c cbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
set(all_link_libs_s sbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
set(all_link_libs_cpp butterflypack_cpp sbutterflypack dbutterflypack cbutterflypack zbutterflypack ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
else()
set(all_link_libs_z zbutterflypack mpi_dummy ${ZFP_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB}  m)
set(all_link_libs_d dbutterflypack mpi_dummy ${ZFP_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
set(all_link_libs_c cbutterflypack mpi_dummy ${ZFP_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
set(all_link_libs_s sbutterflypack mpi_dummy ${ZFP_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
set(all_link_libs_cpp butterflypack_cpp sbutterflypack dbutterflypack cbutterflypack zbutterflypack mpi_dummy ${ZFP_LIB} ${MAGMA_LIB} ${LAPACK_LIB} ${BLAS_LIB} m)
endif()
########################################### single complex Fortran drivers

set(FCEXM1 EMSURF_Driver_sp.f90 EMSURF_Module_sp.f90)
add_executable(ie3d_sp ${FCEXM1})
target_link_libraries(ie3d_sp ${all_link_libs_c})

########################################### double complex Fortran drivers

set(FZEXM1 EMCURV_Driver.f90 EMCURV_Module.f90)
add_executable(ie2d ${FZEXM1})
target_link_libraries(ie2d ${all_link_libs_z})
install(TARGETS ie2d RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")
add_dependencies(ie2d ie3d_sp)

set(FZEXM2 EMCURV_Eigen_Driver.f90 EMCURV_Module.f90)
add_executable(ie2deigen ${FZEXM2})
target_link_libraries(ie2deigen ${all_link_libs_z})
add_dependencies(ie2deigen ie2d)

set(FZEXM3 EMSURF_Driver.f90 EMSURF_Module.f90)
add_executable(ie3d ${FZEXM3})
target_link_libraries(ie3d ${all_link_libs_z})
add_dependencies(ie3d ie2deigen)

set(FZEXM4 EMSURF_Eigen_Driver.f90 EMSURF_Module.f90)
add_executable(ie3deigen ${FZEXM4})
target_link_libraries(ie3deigen ${all_link_libs_z})
add_dependencies(ie3deigen ie3d)
install(TARGETS ie3deigen RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")

set(FZEXM5 SMAT_Driver.f90)
add_executable(smat ${FZEXM5})
target_link_libraries(smat ${all_link_libs_z})
add_dependencies(smat ie3deigen)

set(FZEXM6 FrontalDist_Driver.f90)
add_executable(frontaldist ${FZEXM6})
target_link_libraries(frontaldist ${all_link_libs_z})
add_dependencies(frontaldist smat)

set(FZEXM9 EMSURF_Port_Eigen_Driver.f90 EMSURF_Port_Module.f90 DoCubicInterp2D.c)
add_executable(ie3dporteigen ${FZEXM9})
target_link_libraries(ie3dporteigen ${all_link_libs_z})
add_dependencies(ie3dporteigen frontaldist)
install(TARGETS ie3dporteigen RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")


set(FZEXM10 EMSURF_Port_Driver.f90 EMSURF_Port_Module.f90 DoCubicInterp2D.c)
add_executable(ie3dport ${FZEXM10})
target_link_libraries(ie3dport ${all_link_libs_z})
add_dependencies(ie3dport ie3dporteigen)
install(TARGETS ie3dport RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")


set(FZEXM11 RankBenchmark_Driver.f90)
add_executable(frankben ${FZEXM11})
target_link_libraries(frankben ${all_link_libs_z})
add_dependencies(frankben ie3dport)
install(TARGETS frankben RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")


set(FZEXM12 RankBenchmark_tensor_Driver.f90)
add_executable(frankben_t ${FZEXM12})
target_link_libraries(frankben_t ${all_link_libs_z})
add_dependencies(frankben_t frankben)
install(TARGETS frankben_t RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")

set(FZEXM13 Test_TT_Driver.f90)
add_executable(testtt ${FZEXM13})
target_link_libraries(testtt ${all_link_libs_z})
add_dependencies(testtt frankben_t)
install(TARGETS testtt RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")


set(FZEXM14 Test_QTT_Driver.f90)
add_executable(testqtt ${FZEXM14})
target_link_libraries(testqtt ${all_link_libs_z})
add_dependencies(testqtt testtt)
install(TARGETS testqtt RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")

set(FZEXM15 Test_proxy_selection.f90)
add_executable(testproxy ${FZEXM15})
target_link_libraries(testproxy ${all_link_libs_z})
add_dependencies(testproxy testqtt)
install(TARGETS testproxy RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/EXAMPLE")





########################################### double Fortran drivers
set(FDEXMSEQ1 KERREG_Driver_seq.f90)
add_executable(krr_seq ${FDEXMSEQ1})
target_link_libraries(krr_seq ${all_link_libs_d})

set(FDEXM1 FULLMAT_Driver.f90)
add_executable(full ${FDEXM1})
target_link_libraries(full ${all_link_libs_d})
target_link_libraries(full ${all_link_libs_z})
add_dependencies(full krr_seq)

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


########################################### non-MPI double C++ drivers
set(CDEXM_SEQ InterfaceTest_simple_seq.cpp)
add_executable(ctest_simple_seq ${CDEXM_SEQ})
target_link_libraries(ctest_simple_seq ${all_link_libs_d})

set(CDEXM_SEQ1 InterfaceTest_simple_seq_newapi.cpp)
add_executable(ctest_simple_seq_newapi ${CDEXM_SEQ1})
target_link_libraries(ctest_simple_seq_newapi ${all_link_libs_cpp})

if(enable_mpi)

########################################### double C++ drivers
set(CDEXM1 InterfaceTest.cpp)
add_executable(ctest ${CDEXM1})
target_link_libraries(ctest ${all_link_libs_d})

set(CDEXM2 InterfaceTest_simple.cpp)
add_executable(ctest_simple ${CDEXM2})
target_link_libraries(ctest_simple ${all_link_libs_d})
add_dependencies(ctest_simple ctest)

set(CDEXM3 InterfaceTest_simple_newapi.cpp)
add_executable(ctest_simple_newapi ${CDEXM3})
target_link_libraries(ctest_simple_newapi ${all_link_libs_cpp})


########################################### double complex C++ drivers

set(CZEXM1 Taylor2D.cpp G2D/G2D.cpp G2D/rtGetInf.cpp G2D/rtGetNaN.cpp G2D/rt_nonfinite.cpp G2D/bessel.cpp)
add_executable(go2d ${CZEXM1})
target_link_libraries(go2d ${all_link_libs_z})
add_dependencies(go2d ctest_simple)

set(CZEXM2 Taylor3D.cpp G3D/G3D.cpp G3D/rtGetInf.cpp G3D/rtGetNaN.cpp G3D/rt_nonfinite.cpp)
add_executable(go3d ${CZEXM2})
target_link_libraries(go3d ${all_link_libs_z})
add_dependencies(go3d go2d)

# exclude Clang as Clang<13 doesn't support creal
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CZEXM3 FIO_Driver.cpp)
add_executable(cfio ${CZEXM3})
target_link_libraries(cfio ${all_link_libs_d})
target_link_libraries(cfio ${all_link_libs_z})
add_dependencies(cfio go3d)

set(CZEXM4 RankBenchmark_Driver.cpp)
add_executable(rankben ${CZEXM4})
target_link_libraries(rankben ${all_link_libs_d})
target_link_libraries(rankben ${all_link_libs_z})
add_dependencies(rankben cfio)
endif()

# exclude Clang as Clang<13 doesn't support cabs,conj
if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CIFEXM InverseFIO_Driver.cpp)
add_executable(cifio ${CIFEXM})
target_link_libraries(cifio ${all_link_libs_d})
target_link_libraries(cifio ${all_link_libs_z})
add_dependencies(cifio rankben)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CIFEXM2D InverseFIO2D_Driver.cpp)
add_executable(cifio2d ${CIFEXM2D})
target_link_libraries(cifio2d ${all_link_libs_d})
target_link_libraries(cifio2d ${all_link_libs_z})
add_dependencies(cifio2d cifio)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CIFEXM2D InverseFIO2D_SB_Driver.cpp)
add_executable(cifio2dsb ${CIFEXM2D})
target_link_libraries(cifio2dsb ${all_link_libs_d})
target_link_libraries(cifio2dsb ${all_link_libs_z})
add_dependencies(cifio2dsb cifio2d)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CVIE2D VIE2D_Driver.cpp G2D/bessel.cpp)
add_executable(cvie2d ${CVIE2D})
target_link_libraries(cvie2d ${all_link_libs_d})
target_link_libraries(cvie2d ${all_link_libs_z})
add_dependencies(cvie2d cifio2dsb)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CVIE2DT VIE2D_tensor_Driver.cpp G2D/bessel.cpp)
add_executable(cvie2d_t ${CVIE2DT})
target_link_libraries(cvie2d_t ${all_link_libs_d})
target_link_libraries(cvie2d_t ${all_link_libs_z})
add_dependencies(cvie2d_t cvie2d)
endif()


if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CVIE3D VIE3D_Driver.cpp G2D/bessel.cpp)
add_executable(cvie3d ${CVIE3D})
target_link_libraries(cvie3d ${all_link_libs_d})
target_link_libraries(cvie3d ${all_link_libs_z})
add_dependencies(cvie3d cvie2d_t)
endif()


if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
set(CVIE3DT VIE3D_tensor_Driver.cpp G2D/bessel.cpp)
add_executable(cvie3d_t ${CVIE3DT})
target_link_libraries(cvie3d_t ${all_link_libs_d})
target_link_libraries(cvie3d_t ${all_link_libs_z})
add_dependencies(cvie3d_t cvie3d)
endif()


endif()





