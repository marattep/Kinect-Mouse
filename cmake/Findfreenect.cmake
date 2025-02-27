find_path(FREENECT_INCLUDE_DIR libfreenect.h
  PATHS /usr/include /usr/local/include
)

find_library(FREENECT_LIBRARY NAMES freenect
  PATHS /usr/lib /usr/local/lib  /usr/lib/x86_64-linux-gnu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(freenect DEFAULT_MSG FREENECT_LIBRARY FREENECT_INCLUDE_DIR)

if(FREENECT_FOUND)
  set(FREENECT_LIBRARIES ${FREENECT_LIBRARY})
  set(FREENECT_INCLUDE_DIRS ${FREENECT_INCLUDE_DIR})
endif()
