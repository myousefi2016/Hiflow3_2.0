# - Try to find PETSC
#

find_path (PETSC_DIR include/petsc.h HINTS ENV PETSC_DIR DOC "PETSC Directory")

IF(EXISTS ${PETSC_DIR}/include/petsc.h)
  SET(PETSC_FOUND YES)
  SET(PETSC_INCLUDES ${PETSC_DIR})
  find_path (PETSC_INCLUDE_DIR petsc.h HINTS "${PETSC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND PETSC_INCLUDES ${PETSC_INCLUDE_DIR})
  find_library(PETSC_LIBRARY NAMES libpetsc.so PATHS ${PETSC_DIR}/lib)
  list(APPEND PETSC_LIBRARIES ${PETSC_LIBRARY})
ELSE(EXISTS ${PETSC_DIR}/include/petsc.h)
  SET(PETSC_FOUND NO)
ENDIF(EXISTS ${PETSC_DIR}/include/petsc.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDES)
