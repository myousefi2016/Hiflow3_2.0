# - Try to find ILU++
# Once done, this will define
#
#  ILUPP_FOUND - system has ILU++
#  ILUPP_INCLUDE_DIRS - the ILU++ include directories
#  ILUPP_LIBRARY - library to link to

include(LibFindMacros)

# Dependencies
# 

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(ILUPP_PKGCONF ILUPP)

# Include dir
find_path(ILUPP_INCLUDE_DIR
  NAMES iluplusplus_interface.h
  PATHS ${ILUPP_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES iluplusplus
)

# Finally the library itself
# TODO: does the version number have to be included?
find_library(ILUPP_LIBRARY
  NAMES iluplusplus-1.1
  PATHS ${ILUPP_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(ILUPP_PROCESS_INCLUDES ILUPP_INCLUDE_DIR)
set(ILUPP_PROCESS_LIBS ILUPP_LIBRARY)
libfind_process(ILUPP)
