# - Config file for the ButterflyPACK package
#
# It defines the following variables
#  ButterflyPACK_FOUND        - ButterflyPACK was found or not
#  ButterflyPACK_INCLUDE_DIRS - include directories for ButterflyPACK
#  ButterflyPACK_LIBRARIES    - libraries to link against
#
# And the following imported targets:
#   ButterflyPACK::dbutterflypack ButterflyPACK::zbutterflypack
#
include("${CMAKE_CURRENT_LIST_DIR}/butterflypack-config-version.cmake")

include(FindPackageHandleStandardArgs)
set(${CMAKE_FIND_PACKAGE_NAME}_CONFIG "${CMAKE_CURRENT_LIST_FILE}")
find_package_handle_standard_args(${CMAKE_FIND_PACKAGE_NAME} CONFIG_MODE)

if(NOT TARGET ButterflyPACK AND NOT ButterflyPACK_BINARY_DIR)
  include(${CMAKE_CURRENT_LIST_DIR}/butterflypack-targets.cmake)
endif()

set(ButterflyPACK_LIBRARIES ButterflyPACK::dbutterflypack ButterflyPACK::zbutterflypack)
set(ButterflyPACK_INCLUDE_DIRS
  $<TARGET_PROPERTY:ButterflyPACK::dbutterflypack,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:ButterflyPACK::zbutterflypack,INTERFACE_INCLUDE_DIRECTORIES>
)
