# - Try to find BLAS
# Once done, this will define
#
#  BLAS_FOUND - system has BLAS
#  BLAS_INCLUDE_DIRS - the BLAS include directories
#  BLAS_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
#

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(BLAS_PKGCONF BLAS)

# Finally the library itself
find_library(BLAS_LIBRARY
  NAMES libblas.a
  PATHS  ${BLAS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(BLAS_PROCESS_LIBS BLAS_LIBRARY)
libfind_process(BLAS)
