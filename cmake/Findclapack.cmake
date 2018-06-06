# - Try to find CLAPACK
# Once done, this will define
#
#  CLAPACK_FOUND - system has CBLAS
#  CLAPACK_INCLUDE_DIR - the CBLAS include directories
#  CLAPACK_LIBRARIES - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths

# Include dir
find_path(CLAPACK_INCLUDE_DIR
  NAMES lapacke.h
  #PATHS  ${CBLAS_PKGCONF_LIBRARY_DIRS}
)

# Finally the library itself
find_library(CLAPACK_LIBRARIES
  NAMES lapacke
  #PATHS  ${CBLAS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(CLAPACK_PROCESS_INCLUDES CLAPACK_INCLUDE_DIR)
set(CLAPACK_PROCESS_LIBS CLAPACK_LIBRARIES)
libfind_process(CLAPACK)
