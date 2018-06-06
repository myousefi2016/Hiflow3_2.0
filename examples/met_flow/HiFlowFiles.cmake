set(SUBDIRS_LIBRARY
  assembler
  general
  tmp_config
)

# NOTE: Can only compile one application at once, since applications use different config parameters stored in met_flow_vars.h
set(SUBDIRS_EXECUTABLE
    conv_diff
    #wavetank
)

set(METFLOW_BASE_DIR ${PROJECT_SOURCE_DIR}/examples/met_flow)
set(METFLOW_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/met_flow/include)
set(METFLOW_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/met_flow)

# Recurse through sub-directories and gather information for libhiflow
foreach(dir ${SUBDIRS_EXECUTABLE})
  file(GLOB CONFIG_FILES "${METFLOW_BASE_DIR}/${dir}/met_flow_vars.h")
  file(COPY ${CONFIG_FILES} DESTINATION ${METFLOW_BASE_DIR}/tmp_config)
endforeach()

foreach(dir ${SUBDIRS_LIBRARY})
  include(${METFLOW_BASE_DIR}/tmp_config/HiFlowFiles.cmake)

  include(${METFLOW_BASE_DIR}/${dir}/HiFlowFiles.cmake)

  foreach(i ${${dir}_SOURCES})
    list(APPEND METFLOW_SOURCES ${METFLOW_BASE_DIR}/${dir}/${i})
  endforeach()

  foreach(i ${${dir}_PUBLIC_HEADERS})
    # Build directory include file from public headers
    list(APPEND ${dir}_HEADER_FILES "#include \"${dir}/${i}\"")

    # Copy public headers to build directory
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/met_flow/${dir}/${i}" "${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/${dir}/${i}" COPYONLY)
  endforeach()

  string(TOUPPER "${dir}" DIRECTORY)

  # Create directory header file template
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/met_flow/directory.h.in" "${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/${dir}.h.in" @ONLY)

  # Create directory header file
  string(REPLACE ";" "\n" ${DIRECTORY}_HEADER_FILES "${${dir}_HEADER_FILES}")
  configure_file("${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/${dir}.h.in" "${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/${dir}.h")
  list(APPEND METFLOW_HEADER_FILES "#include \"${dir}.h\"")
endforeach()

# Include config.h
list(APPEND METFLOW_HEADER_FILES "#include \"config.h\"")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/")

# Create overall metflow.h include file
string(REPLACE ";" "\n" METFLOW_HEADER_FILES "${METFLOW_HEADER_FILES}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/met_flow/metflow.h.in" "${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/metflow.h")

# External includes
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})

# Build MetFlow library
add_library(metflow ${METFLOW_SOURCES})
set_target_properties(metflow PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${METFLOW_LIB_DIR})

# Add dependencies
target_link_libraries(metflow hiflow ${MPI_LIBRARIES})

# Create Directories
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/in)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/out)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/mesh)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/log)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/start)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/snapshots)

make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/eigvec)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/eigval)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/eigfunc) 
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/obj)

make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/indicators)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/out_patch)
make_directory(${PROJECT_BINARY_DIR}/examples/met_flow/bin/out_inter)

# Make executables
foreach(dir ${SUBDIRS_EXECUTABLE})
  include(${METFLOW_BASE_DIR}/${dir}/HiFlowFiles.cmake)

  add_executable(${dir} ${METFLOW_BASE_DIR}/${dir}/main.cc)
  foreach(i ${${dir}_SOURCES})
    target_sources(${dir} PUBLIC ${METFLOW_BASE_DIR}/${dir}/${i})
  endforeach()

  set_target_properties(${dir} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples/met_flow/bin)
  target_link_libraries(${dir} metflow)
  if (WITH_HDF5)
     target_link_libraries(${dir} hiflow hdf5 dl)
  endif()
  # TODO copy to install directory

  file(GLOB INPUT_FILES "${METFLOW_BASE_DIR}/${dir}/in/*")

  file(COPY ${INPUT_FILES}
    DESTINATION ${PROJECT_BINARY_DIR}/examples/met_flow/bin/in)
endforeach()

# TODO link against METIS, PARMETIS, P4EST, HYPRE, PETSC, SLEPC ?

#TODO
# install library
# install(TARGETS metflow EXPORT metflow
#  ARCHIVE DESTINATION ${METFLOW_LIB_DIR})

# install public header files
#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/met_flow/include/
#  DESTINATION ${HIFLOW_INCLUDE_DIR}
#  FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
#  FILES_MATCHING PATTERN "*.h")
