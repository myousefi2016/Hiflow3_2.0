include_directories(${Boost_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR}/src/include)
include_directories(${MPI_INCLUDE_PATH})

if (WITH_ILUPP)
  include_directories(${ILUPP_INCLUDE_DIR})
endif()

if (WITH_HDF5)
  include_directories(${HDF5_INCLUDE_DIR})
endif()

if (WITH_HYPRE)
  include_directories(${HYPRE_INCLUDE_DIR})
  link_directories(${HYPRE_LIBRARIES})
endif()

# MUMPS direct solver
if (WITH_MUMPS)
  include_directories(${MUMPS_INCLUDE_DIR})
endif()

if (WITH_UMFPACK)
  include_directories(${UMFPACK_INCLUDE_DIR})
endif()

if (WITH_CUDA)
  include_directories(${CUDA_TOOLKIT_INCLUDE})
endif()

if(WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
endif()

if(WITH_P4EST)
  include_directories(${P4EST_INCLUDE_DIR})
endif()

add_executable(inp2vtu inp2vtu.cc)
target_link_libraries(inp2vtu hiflow)
if (WITH_HDF5)
  target_link_libraries(inp2vtu hiflow hdf5 dl)
endif()

add_executable(check_mesh_input check_mesh_input.cc)
target_link_libraries(check_mesh_input hiflow)
if (WITH_HDF5)
  target_link_libraries(check_mesh_input hiflow hdf5 dl)
endif()

add_executable(check_partitioning check_partitioning.cc)
target_link_libraries(check_partitioning hiflow)
if (WITH_HDF5)
  target_link_libraries(check_partitioning hiflow hdf5 dl)
endif()

add_executable(mesh_refiner mesh_refiner.cc)
target_link_libraries(mesh_refiner hiflow)
if (WITH_HDF5)
  target_link_libraries(mesh_refiner hiflow hdf5 dl)
endif()

#add_executable(griddata2vtk griddata2vtk.cc)
#target_link_libraries(griddata2vtk hiflow)

add_executable(numconn number_cell_connections.cc)
target_link_libraries(numconn hiflow)
if (WITH_HDF5)
  target_link_libraries(numconn hiflow hdf5 dl)
endif()

add_executable(quadrature_solin_generator quadrature_solin_generator.cc)

add_executable(partition_mesh partition_mesh.cc)
target_link_libraries(partition_mesh hiflow)
if (WITH_HDF5)
  target_link_libraries(partition_mesh hiflow hdf5 dl)
endif()

install(TARGETS inp2vtu check_mesh_input check_partitioning
  RUNTIME DESTINATION ${HIFLOW_BIN_DIR})

install(TARGETS check_partitioning
  DESTINATION ${HIFLOW_BIN_DIR})

# list of scripts to install
set(UTILS_SCRIPTS info2html.pl mesh2vtu.py pvtu2vtu.py)
# REMOVE BEFORE RELEASE BEGIN
set(UTILS_SCRIPTS info2html.pl gprof2dot.py mesh2vtu.py pvtu2vtu.py)
# REMOVE BEFORE RELEASE END

install(PROGRAMS ${UTILS_SCRIPTS}
  DESTINATION ${HIFLOW_BIN_DIR})
