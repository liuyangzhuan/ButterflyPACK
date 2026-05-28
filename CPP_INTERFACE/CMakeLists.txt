include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_DOUBLECOMPLEX)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_SINGLE)
include_directories(${ButterflyPACK_SOURCE_DIR}/SRC_COMPLEX)
include_directories(${ButterflyPACK_SOURCE_DIR}/EXAMPLE)
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
    m_BPACK_utilities.f90
)



add_library(butterflypack_cpp OBJECT ${sources} ${HEADERS})
add_dependencies(butterflypack_cpp
  dbutterflypack
  zbutterflypack
  sbutterflypack
  cbutterflypack)

set_target_properties(butterflypack_cpp PROPERTIES
  POSITION_INDEPENDENT_CODE ON)

target_include_directories(butterflypack_cpp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

# Define GNU standard installation directories
include(GNUInstallDirs)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


#SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR})
# install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
