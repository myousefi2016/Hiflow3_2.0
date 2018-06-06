# - Try to find UMFPACK
# Once done, this will define
#
#  UMFPACK_FOUND - system has UMFPACK
#  UMFPACK_INCLUDE_DIRS - the UMFPACK include directories
#  UMFPACK_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(UMFPACK_PKGCONF UMFPACK)

# Include dir
find_path(UMFPACK_INCLUDE_DIR
  NAMES umfpack.h
  PATHS ${UMFPACK_PKGCONF_INCLUDE_DIRS}
)

find_library(AMD_LIBRARY
  NAMES libamd.a
  PATHS ${UMFPACK_PKGCONF_LIBRARY_DIRS}
)

# Finally the library itself
find_library(UMFPACK_LIBRARY
  NAMES libumfpack.a
  PATHS ${UMFPACK_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(UMFPACK_PROCESS_INCLUDES UMFPACK_INCLUDE_DIR)
set(UMFPACK_PROCESS_LIBS "UMFPACK_LIBRARY")
libfind_process(UMFPACK)
