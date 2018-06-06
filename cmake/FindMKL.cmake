# - Try to find MKL
# Once done, this will define
#
#  MKL_FOUND - system has MKL
#  MKL_INCLUDE_DIRS - the MKL include directories
#  MKL_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(MKL_PKGCONF MKL)

# Include dir
find_path(MKL_INCLUDE_DIR
  NAMES mkl.h
)

# Finally the library itself
find_library(MKL_LIBRARY
  NAMES libmkl_intel_lp64.a
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(MKL_PROCESS_INCLUDES MKL_INCLUDE_DIR)
set(MKL_PROCESS_LIBS MKL_LIBRARY)
libfind_process(MKL)
