include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(SUBDIRS
  assembly
  common
  dof
  eigen_value
  fem
  linear_algebra
  linear_solver
  mesh
  nonlinear
  polynomial_chaos
  polynomials
  quadrature
  space
  tools
  adaptivity
)

# Recurse through sub-directories and gather information for libhiflow
foreach(dir ${SUBDIRS})

  include(${dir}/HiFlowFiles.cmake)

  foreach(i ${${dir}_SOURCES})
    list(APPEND HIFLOW_SOURCES ${dir}/${i})
  endforeach()

  foreach(i ${${dir}_PUBLIC_HEADERS})
    # Build directory include file from public headers
    list(APPEND ${dir}_HEADER_FILES "#include \"${dir}/${i}\"")

    # Copy public headers to build directory
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${i}" "${CMAKE_CURRENT_BINARY_DIR}/include/${dir}/${i}" COPYONLY)
  endforeach()

  string(TOUPPER "${dir}" DIRECTORY)

  # Create directory header file template
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/directory.h.in" "${CMAKE_CURRENT_BINARY_DIR}/include/${dir}.h.in" @ONLY)

  # Create directory header file
  string(REPLACE ";" "\n" ${DIRECTORY}_HEADER_FILES "${${dir}_HEADER_FILES}")
  configure_file("${CMAKE_CURRENT_BINARY_DIR}/include/${dir}.h.in" "${CMAKE_CURRENT_BINARY_DIR}/include/${dir}.h")
  list(APPEND HIFLOW_HEADER_FILES "#include \"${dir}.h\"")

endforeach()

# Include config.h
list(APPEND HIFLOW_HEADER_FILES "#include \"config.h\"")


include_directories("${CMAKE_CURRENT_BINARY_DIR}/include/")

# Create overall hiflow.h include file
string(REPLACE ";" "\n" HIFLOW_HEADER_FILES "${HIFLOW_HEADER_FILES}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/hiflow.h.in" "${CMAKE_CURRENT_BINARY_DIR}/include/hiflow.h")

# External includes
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR}/contrib/tinyxml)

# Add tinyxml source to be compiled in library.
set(TINYXML_DIR ${PROJECT_SOURCE_DIR}/contrib/tinyxml)
set(TINYXML_SOURCES tinyxml.cpp tinyxmlparser.cpp tinyxmlerror.cpp tinystr.cpp)
foreach(i ${TINYXML_SOURCES})
  list(APPEND HIFLOW_SOURCES "${TINYXML_DIR}/${i}")
endforeach()

# Compile CUDA files
if (WITH_CUDA)
  include_directories(${CUDA_TOOLKIT_INCLUDE})
  foreach(i ${HIFLOW_SOURCES})

    if (${i} MATCHES ".*\\.cu$")
      cuda_compile(c ${i})
      list(REMOVE_ITEM HIFLOW_SOURCES ${i})
      list(APPEND HIFLOW_SOURCES ${c})
    endif()

  endforeach()
endif()

if (WITH_ILUPP)
  include_directories(${ILUPP_INCLUDE_DIR})
endif()

if (WITH_UMFPACK)
  include_directories(${UMFPACK_INCLUDE_DIR})
endif()

# Build HiFlow library
add_library(hiflow ${HIFLOW_SOURCES})

# Add dependencies
target_link_libraries(hiflow ${MPI_LIBRARIES})

# OpenCL
if (WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
  target_link_libraries(hiflow ${OPENCL_LIBRARIES})
endif()


# HDF5
if (WITH_HDF5)
  target_link_libraries(hiflow hdf5)
  include_directories(${HDF5_INCLUDE_DIR})
  link_directories(${HDF5_LIBRARY_DIRS})
endif()

# MKL BLAS/LAPACK
if (WITH_MKL)
  target_link_libraries(hiflow mkl)
  include_directories(${MKL_INCLUDE_DIR})
  #link_directories(${MKL_LIBRARY})
endif()

# OTHER CBLAS library
if (WITH_CBLAS)
  include_directories(${CBLAS_INCLUDE_DIR})
  target_link_libraries(hiflow ${CBLAS_LIBRARY})
endif()

if (WITH_CLAPACK)
  include_directories(${CLAPACK_INCLUDE_DIR})
  target_link_libraries(hiflow ${CLAPACK_LIBRARIES})
endif()

# MUMPS direct solver
if (WITH_MUMPS)
  include_directories(${MUMPS_INCLUDE_DIR})
  target_link_libraries(hiflow mumps)
  target_link_libraries(hiflow pord)
  target_link_libraries(hiflow scalapack)
  target_link_libraries(hiflow blacs)
  target_link_libraries(hiflow blacs_init)
  target_link_libraries(hiflow blas)
  target_link_libraries(hiflow metis)
endif()

# UMFPACK sequential direct solver
if (WITH_UMFPACK)
  target_link_libraries(hiflow umfpack)
endif()

# CUDA libraries
if (WITH_CUDA)
  target_link_libraries(hiflow ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY})
endif()

# GAUSSQ Fortran library
if (WITH_GAUSSQ)
  target_link_libraries(hiflow ${PROJECT_SOURCE_DIR}/contrib/d1mach.f ${PROJECT_SOURCE_DIR}/contrib/gaussq.f -lgfortran)
endif()

# ILU++ library
if (WITH_ILUPP)
  include_directories(${ILUPP_INCLUDE_DIR})
  target_link_libraries(hiflow ilupp)
endif()

# METIS graph partitioning
# Must come after ILU++ due to dependencies
if (WITH_METIS)
  include_directories(${METIS_INCLUDE_DIR})
  target_link_libraries(hiflow metis)
endif()

# ParMETIS partitioning
if (WITH_PARMETIS)
  include_directories(${ParMETIS_INCLUDE_DIR})
  target_link_libraries(hiflow parmetis)
endif()

# HYPRE interface
if(WITH_HYPRE)
  include_directories(${HYPRE_INCLUDE_DIR})
  target_link_libraries(hiflow HYPRE)
endif()

# PETSC 
if(WITH_PETSC)
  include_directories(${PETSC_INCLUDE_DIR})
  target_link_libraries(hiflow PETSC)
endif()
if(WITH_COMPLEX_PETSC)
  include_directories(${PETSC_INCLUDE_DIR})
  target_link_libraries(hiflow PETSC)
endif()

# SLEPC 
if(WITH_SLEPC)
  include_directories(${SLEPC_INCLUDE_DIR})
  target_link_libraries(hiflow SLEPC)
endif()

# P4EST
if(WITH_P4EST)
  include_directories(${P4EST_INCLUDE_DIR})
  target_link_libraries(hiflow P4EST)
  target_link_libraries(hiflow P4EST_SC)
  target_link_libraries(hiflow P4EST_Z)
endif()

# install library
install(TARGETS hiflow EXPORT HiFlow
  ARCHIVE DESTINATION ${HIFLOW_LIB_DIR})

# install public header files
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/
  DESTINATION ${HIFLOW_INCLUDE_DIR}
  FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
  FILES_MATCHING PATTERN "*.h")
