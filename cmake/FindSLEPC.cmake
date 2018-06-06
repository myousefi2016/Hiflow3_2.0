# - Try to find SLEPC
#

find_path (SLEPC_DIR include/slepc.h HINTS ENV SLEPC_DIR DOC "SLEPC Directory")

IF(EXISTS ${SLEPC_DIR}/include/slepc.h)
  SET(SLEPC_FOUND YES)
  SET(SLEPC_INCLUDES ${SLEPC_DIR})
  find_path (SLEPC_INCLUDE_DIR slepc.h HINTS "${SLEPC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND SLEPC_INCLUDES ${SLEPC_INCLUDE_DIR})
  find_library(SLEPC_LIBRARY NAMES libslepc.so PATHS ${SLEPC_DIR}/lib)
  list(APPEND SLEPC_LIBRARIES ${SLEPC_LIBRARY})
ELSE(EXISTS ${SLEPC_DIR}/include/slepc.h)
  SET(SLEPC_FOUND NO)
ENDIF(EXISTS ${SLEPC_DIR}/include/slepc.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPC DEFAULT_MSG SLEPC_LIBRARIES SLEPC_INCLUDES)
