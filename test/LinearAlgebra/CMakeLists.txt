set(TestName LinearAlgebra)

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/tools
    ${PROJECT_BINARY_DIR}/src/include
    ${Boost_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
    ${GTEST_INCLUDE_DIR}
)

if (WITH_ILUPP)
  include_directories(${ILUPP_INCLUDE_DIR})
endif()

if (WITH_MUMPS)
  include_directories(${MUMPS_INCLUDE_DIR})
endif()

if (WITH_UMFPACK)
  include_directories(${UMFPACK_INCLUDE_DIR})
endif()

if(WITH_CBLAS)
  include_directories(${CBLAS_INCLUDE_DIR})
  link_directories(${CBLAS_LIBRARY})
endif()

if(WITH_CLAPACK)
  include_directories(${CBLAS_INCLUDE_DIR})
  include_directories(${CLAPACK_INCLUDE_DIR})
  link_directories(${CLAPACK_LIBRARIES})
  link_directories(${CBLAS_LIBRARY})
endif()

if (WITH_MKL)
  include_directories(${MKL_INCLUDE_DIR})
  link_directories(${MKL_LIBRARY})
endif()

if (WITH_HDF5)
  include_directories(${HDF5_INCLUDE_DIR})
endif()

if (WITH_HYPRE)
  include_directories(${HYPRE_INCLUDE_DIR})
  link_directories(${HYPRE_LIBRARIES})
endif()

if (WITH_PETSC)
  include_directories(${PETSC_INCLUDE_DIR})
endif()

if (WITH_SLEPC)
  include_directories(${SLEPC_INCLUDE_DIR})
endif()

if(WITH_COMPLEX_PETSC)
  include_directories(${PETSC_INCLUDE_DIR})
endif()

if(WITH_P4EST)
  include_directories(${P4EST_INCLUDE_DIR})
endif()

if(WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
endif()

add_executable(
    ${TestName}
    cg.cc
    gmres.cc
    main.cc
    #solver.cc
    CPUsimple_CSR_lMatrix.cc
)

target_link_libraries(
    ${TestName}
    hiflow
    ${GTEST_BOTH_LIBRARIES}
    -pthread
)

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 1 ${PROJECT_BINARY_DIR}/test/LinearAlgebra/${TestName} --gtest_color=no --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)
