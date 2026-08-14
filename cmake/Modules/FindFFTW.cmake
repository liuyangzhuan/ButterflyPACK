# cmake/FindFFTW.cmake

find_path(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  HINTS $ENV{FFTW_ROOT}
  PATH_SUFFIXES include
)

find_library(FFTW_DOUBLE_LIBRARY
  NAMES fftw3
  HINTS $ENV{FFTW_ROOT}
  PATH_SUFFIXES lib lib64
)

find_library(FFTW_SINGLE_LIBRARY
  NAMES fftw3f
  HINTS $ENV{FFTW_ROOT}
  PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW
  REQUIRED_VARS
    FFTW_INCLUDE_DIR
    FFTW_DOUBLE_LIBRARY
    FFTW_SINGLE_LIBRARY
)

if(FFTW_FOUND)
  set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
  set(FFTW_LIBRARIES ${FFTW_DOUBLE_LIBRARY} ${FFTW_SINGLE_LIBRARY})

  if(NOT TARGET FFTW::Double)
    add_library(FFTW::Double UNKNOWN IMPORTED)
    set_target_properties(FFTW::Double PROPERTIES
      IMPORTED_LOCATION "${FFTW_DOUBLE_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
    )
  endif()

  if(NOT TARGET FFTW::Single)
    add_library(FFTW::Single UNKNOWN IMPORTED)
    set_target_properties(FFTW::Single PROPERTIES
      IMPORTED_LOCATION "${FFTW_SINGLE_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
    )
  endif()

  if(NOT TARGET FFTW::FFTW)
    add_library(FFTW::FFTW INTERFACE IMPORTED)
    target_link_libraries(FFTW::FFTW INTERFACE FFTW::Double FFTW::Single)
  endif()
endif()

mark_as_advanced(
  FFTW_INCLUDE_DIR
  FFTW_DOUBLE_LIBRARY
  FFTW_SINGLE_LIBRARY
)