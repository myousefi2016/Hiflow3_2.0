# - Try to find ParMETIS
# Once done, this will define
#
#  ParMETIS_FOUND - system has ParMETIS
#  ParMETIS_INCLUDE_DIRS - the ParMETIS include directories
#  ParMETIS_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(ParMETIS_PKGCONF ParMETIS)

# Include dir
find_path(ParMETIS_INCLUDE_DIR
  NAMES parmetis.h
  PATHS ${ParMETIS_PKGCONF_INCLUDE_DIRS}
)

# Only check for the library since this is a link dependency
find_library(ParMETIS_LIBRARY
  NAMES parmetis
  PATHS ${ParMETIS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(ParMETIS_PROCESS_INCLUDES ParMETIS_INCLUDE_DIR)
set(ParMETIS_PROCESS_LIBS ParMETIS_LIBRARY)
libfind_process(ParMETIS)
