set(TestName AMG)
add_definitions(-DDATA_DIR="${PROJECT_SOURCE_DIR}/test/AMG/data/" -DGTEST_USE_OWN_TR1_TUPLE=0)

include_directories(
    ${GTEST_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/tools
    ${PROJECT_BINARY_DIR}/src/include
    ${Boost_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
)
if (WITH_ILUPP)
    include_directories(${ILUPP_INCLUDE_DIR})
endif()

if(WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
endif()

add_executable(
    ${TestName}
    DirectInterpolation.cc
    intersection_size.cc
    LevelGenerator.cc
    main.cc
    ReturnValueOptimization.cc
    StandardCoarsening.cc
    StandardInterpolation.cc
)

target_link_libraries(
    ${TestName}
    hiflow
    ${GTEST_BOTH_LIBRARIES}
    pthread
)
if (WITH_HDF5)
    target_link_libraries(
        ${TestName}
        hiflow
        ${GTEST_BOTH_LIBRARIES}
        pthread hdf5 dl
    )
endif()

add_test(
    NAME ${TestName}
    COMMAND ${MPIEXEC} -np 1 ${PROJECT_BINARY_DIR}/test/AMG/${TestName} --gtest_color=no --gtest_output=xml:${CMAKE_BINARY_DIR}/Testing/${TestName}.xml
)
