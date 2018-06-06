# - Try to find CBLAS
# Once done, this will define
#
#  CBLAS_FOUND - system has CBLAS
#  CBLAS_INCLUDE_DIRS - the CBLAS include directories
#  CBLAS_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(CBLAS_PKGCONF CBLAS)

# Include dir
find_path(CBLAS_INCLUDE_DIR
  NAMES cblas.h
  PATHS  ${CBLAS_PKGCONF_LIBRARY_DIRS}
)

# Finally the library itself
find_library(CBLAS_LIBRARY
  NAMES libcblas.a libatlas.a
  PATHS  ${CBLAS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(CBLAS_PROCESS_INCLUDES CBLAS_INCLUDE_DIR)
set(CBLAS_PROCESS_LIBS CBLAS_LIBRARY)
libfind_process(CBLAS)
