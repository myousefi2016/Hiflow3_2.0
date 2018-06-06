set(ExecutableName CoupledMatrix_MPI)

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

if(WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
endif()

add_executable(
    ${ExecutableName}
    CoupledMatrix_dense.cc
    CoupledMatrix_mesh.cc
    main.cc
)

target_link_libraries(
    ${ExecutableName}
    hiflow
    ${GTEST_BOTH_LIBRARIES}
    pthread
)
if (WITH_HDF5)
    target_link_libraries(${ExecutableName}
    hiflow
    ${GTEST_BOTH_LIBRARIES}
    pthread hdf5 dl)
  endif()

set(TestName ${ExecutableName}_2cores)

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 2 ${PROJECT_BINARY_DIR}/test/CoupledMatrix_MPI/${ExecutableName} --gtest_color=no --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)

set(TestName ${ExecutableName}_3cores)

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 3 ${PROJECT_BINARY_DIR}/test/CoupledMatrix_MPI/${ExecutableName} --gtest_color=no --gtest_filter=-*CoupledMatrix_mesh.* --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)
