if(NOT enable_mpi)
include_directories(${ButterflyPACK_SOURCE_DIR}/MPI_DUMMY)
endif()
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_DOUBLECOMPLEX)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_SINGLE)
include_directories( ${ButterflyPACK_BINARY_DIR}/SRC_COMPLEX)


set(headers
    CPPWrapper.hpp
)

set(sources
    CPPWrapper.cpp
)



add_library(butterflypack_cpp ${sources} ${HEADERS})

if (enable_mpi)
target_link_libraries(butterflypack_cpp zbutterflypack cbutterflypack dbutterflypack sbutterflypack
                       ${MPI_Fortran_LIBRARIES} ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${SCALAPACK_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
else()
target_link_libraries(butterflypack_cpp mpi_dummy zbutterflypack cbutterflypack dbutterflypack sbutterflypack ${BLAS_LIB} ${LAPACK_LIB} ${ZFP_LIB} ${ARPACK_LIB} ${MAGMA_LIB} m)
endif()
                      set_target_properties(butterflypack_cpp PROPERTIES
                      VERSION ${PROJECT_VERSION} SOVERSION ${VERSION_MAJOR} 
)

target_include_directories(butterflypack_cpp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# Define GNU standard installation directories
include(GNUInstallDirs)

install(TARGETS butterflypack_cpp EXPORT MYBPACKTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


#SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.mod")