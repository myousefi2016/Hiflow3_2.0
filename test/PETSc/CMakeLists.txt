set(ExecutableName PETScTest)
add_definitions(-DDATA_DIR="${PROJECT_SOURCE_DIR}/test/PETSc/data/" -DGTEST_USE_OWN_TR1_TUPLE=0)

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/tools
    ${PROJECT_BINARY_DIR}/src/include
    ${Boost_INCLUDE_DIR}
    ${PETSC_INCLUDE_DIR}
    ${GTEST_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
)

if (WITH_ILUPP)
    include_directories(${ILUPP_INCLUDE_DIR})
endif()

if (WITH_HDF5)
    include_directories(${HDF5_INCLUDE_DIR})
endif()

if(WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
endif()

add_executable(
    ${ExecutableName}
    PETScLinearSolver.cc
    PETScLinearSolverReadData.cc
    PETScMatrix.cc
    PETScVector.cc
    main.cc
)

target_link_libraries(
    ${ExecutableName}
    hiflow
    ${PETSC_LIBRARIES}
    ${GTEST_BOTH_LIBRARIES}
    pthread
)

if (WITH_HDF5)
    target_link_libraries(${ExecutableName}
    hiflow
    ${PETSC_LIBRARIES}
    ${GTEST_BOTH_LIBRARIES}
    pthread hdf5 dl)
endif()

set(TestName ${ExecutableName})

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 1 ${PROJECT_BINARY_DIR}/test/PETSc/${ExecutableName} --gtest_color=no --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)

set(TestName ${ExecutableName}_2cores)

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 2 ${PROJECT_BINARY_DIR}/test/PETSc/${ExecutableName} --gtest_color=no --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)

set(TestName ${ExecutableName}_3cores)

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 3 ${PROJECT_BINARY_DIR}/test/PETSc/${ExecutableName} --gtest_color=no --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)
