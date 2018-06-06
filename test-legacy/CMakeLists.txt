# ======= TESTS ========
#
# This file contains the programs that make up the test suite.  These
# should all verify that some aspect of the various modules function
# correctly.  If a test succeeds, it should return 0, and otherwise it
# should return -1. The include test.h includes some useful C macros
# for writing tests.
#
# The naming convention of tests is that they should exist in a source file with the suffix
# _test.cc . The name of the test is the source filename without this suffix.
#

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${PROJECT_BINARY_DIR}/src/include)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR}/contrib/unittest++)

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

link_directories(${PROJECT_BINARY_DIR}/contrib/unittest++)

# Common defines
add_definitions(-DMESH_DATADIR="${PROJECT_SOURCE_DIR}/test-legacy/data/")
add_definitions(-DHIFLOW_DIR="${PROJECT_SOURCE_DIR}")

# This macros can be used to define a test. It takes as input the name
# of the test (without the _test.cc suffix), and the HiFlow module
# libraries that the test depends on in a space-delimited string. It
# creates a new executable target and adds it to the list of tests
# that are executed when "make test" is run.
macro(def_test TEST_NAME)
  add_executable(${TEST_NAME}_test ${TEST_NAME}_test.cc)
  target_link_libraries(${TEST_NAME}_test hiflow)
  if (WITH_HDF5)
    target_link_libraries(${TEST_NAME}_test hiflow pthread hdf5 dl)
  endif()
  add_test(${TEST_NAME} ${TEST_NAME}_test)
endmacro(def_test)

# TODO: check which of these are real tests, in that they can fail...
# "false" tests are marked with NOTEST

# Common tests
def_test(mpi)                 # NOTEST
def_test(vector_algebra)      # TODO: implement tests

# Dof tests
def_test(dof_identification)
def_test(dof_extended_identification)
def_test(parallel_adaptive_jump_term) # NOTEST

# FEM tests
def_test(fe_weight_ordering)
def_test(fe_weight)
def_test(fe_weight_zero)
def_test(transformation)      # NOTEST

# Mesh tests
def_test(attributes)
def_test(boundary_extraction)
def_test(cell_type)
def_test(incidence_iteration)
def_test(material_number)
def_test(mesh_interface)
def_test(parallel_vtk_reader) # NOTEST
def_test(closest_point)

def_test(reader)              # NOTEST
def_test(vertex_lookup)       # NOTEST
def_test(vtk_reader)          # NOTEST
def_test(writer)              # NOTEST

# Quadrature tests
def_test(quadrature_order)
def_test(facet_integral)

# Vector space tests
def_test(vector_space)
def_test(eval_fe_solution)
def_test(binary_io)
def_test(tensor)

# Assembly tests
add_executable(hp_assembly_test hp_assembly_test.cc)
target_link_libraries(hp_assembly_test hiflow)
if (WITH_HDF5)
  target_link_libraries(hp_assembly_test hiflow hdf5 dl)
endif()
add_test(hp_assembly_test ${PROJECT_BINARY_DIR}/test-legacy/hp_assembly_test 2 1)
def_test(facet_normal) # needs to be fixed

# Shared elements tests
def_test(unittestpp)
target_link_libraries(unittestpp_test UnitTest++)

# TODO needs update
#add_executable(split_mesh_test split_mesh_test.cc)
#target_link_libraries(split_mesh_test hiflow)
#add_test(split_mesh ${MPIEXEC} -np 3 ${PROJECT_BINARY_DIR}/test-legacy/split_mesh_test)

if (WITH_HDF5)
  def_test(hdf5)
endif()

if (RUN_TEDIOUS_TESTS)
  def_test(big_mesh)
  def_test(integration)
  def_test(bdry_integration)
  def_test(refinement)       # FEATURE 90
  def_test(tetra_refinement) # FEATURE 91
endif()

# lmp
def_test(lmp)
target_link_libraries(lmp_test UnitTest++)

# CG
def_test(cg)
target_link_libraries(cg_test UnitTest++)

# CUDA
if(WITH_CUDA)
  def_test(cuda)
  target_link_libraries(cuda_test UnitTest++ ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
endif()

add_executable(point_evaluator_test2D point_evaluator_test2D.cc)
target_link_libraries(point_evaluator_test2D hiflow)
if (WITH_HDF5)
  target_link_libraries(point_evaluator_test2D hiflow hdf5 dl)
endif()
add_test(point_evaluator2D ${MPIEXEC} -np 2 ${PROJECT_BINARY_DIR}/test-legacy/point_evaluator_test2D)

add_executable(point_evaluator_test3D point_evaluator_test3D.cc)
target_link_libraries(point_evaluator_test3D hiflow)
if (WITH_HDF5)
  target_link_libraries(point_evaluator_test3D hiflow hdf5 dl)
endif()
add_test(point_evaluator3D ${MPIEXEC} -np 2 ${PROJECT_BINARY_DIR}/test-legacy/point_evaluator_test3D)

add_executable(incident_vertex_iterator_test incident_vertex_iterator_test.cc)
target_link_libraries(incident_vertex_iterator_test hiflow)
if (WITH_HDF5)
  target_link_libraries(incident_vertex_iterator_test "hiflow" hdf5 dl)
endif()
add_test(incident_vertex_iterator ${MPIEXEC} -np 4 ${PROJECT_BINARY_DIR}/test-legacy/incident_vertex_iterator_test)

#def_demo(memory_test "mesh;dof;fem;space")  # NOTEST
#def_demo(memory_test2 "mesh;dof;fem;space") # NOTEST

# # ===== DEMOS ======
# #
# # Demos are programs that illustrate how different
# # aspects of HiFlow are used.

# # Helper macro to define a demo, whose source code exists in
# # ${DEMO_NAME}.cc, and which depends on the HiFlow modules listed in
# # ${DEMO_LIBS}.
# macro(def_demo DEMO_NAME DEMO_LIBS)
#   add_executable(${DEMO_NAME} ${DEMO_NAME}.cc)
#   foreach(LIB ${DEMO_LIBS})
#     target_link_libraries(${DEMO_NAME} ${LIB})
#   endforeach(LIB)
#   target_link_libraries(${DEMO_NAME} "log")
# endmacro(def_demo)

# # List of demos for Mesh module
# set(MESH_DEMOS
#   iteration_demo
#   advanced_iteration_demo
#   write_refined_mesh
#   cell_type_demo
#   uniform_refinement_perf_test
#   visualization_of_random_data_demo
#   refine_by_attribute
#   local_refinements_by_attribute
#   )

# foreach (DEMO ${MESH_DEMOS})
#   def_demo(${DEMO} "mesh")
# endforeach(DEMO)

# # List of demos for Mesh module that require MPI
# set(MESH_MPI_DEMOS
#   vertex_sync_demo
#   shared_vertex_table_demo
#   ghost_cells_demo
#   )
# # List of demos for dof module that require MPI
# set(DOF_MPI_DEMOS
#   dof_partition_demo
#   )

# # List of demos for Mesh module that require MPI and METIS
# set(MESH_METIS_MPI_DEMOS
#   global_partitioning_demo
#   metis_lung_partitioning_demo
#   )

# include_directories(${MPI_INCLUDE_PATH})
# foreach (DEMO ${MESH_MPI_DEMOS})
#   def_demo(${DEMO} "mesh")
# endforeach (DEMO)

# # special demo requiring mesh tools from applications library
# def_demo(partitioning_demo "utils;mesh")
# def_demo(mesh_tools_partitioning_demo "utils;mesh")
# def_demo(refinement_demo "utils;mesh")
# def_demo(assembly_demo "applications;dof;fem;linear_algebra;mesh;quadrature;space;timer")

# if (WITH_METIS)
#   foreach (DEMO ${MESH_METIS_MPI_DEMOS})
#     def_demo(${DEMO} "mesh")
#   endforeach (DEMO)
# endif (WITH_METIS)

# # Vector space demo
# def_demo(vector_space_demo "mesh;space;fem;dof")
# # dof partition demo
# def_demo(dof_partition_demo "mesh;utils;fem;dof")

# # Laplace demos
# def_demo(laplace2d_convergence_test "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(laplace3d_convergence_test "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(laplace3d_convergence_test_acc "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space;toolbox")
#add_executable(laplace2d_convergence_test_acc laplace2d_convergence_test_acc.cc)
#target_link_libraries(laplace2d_convergence_test_acc hiflow)
# def_demo(convdiff_subalpha_acc "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space;toolbox")
# def_demo(laplace2d_neumann_tutorial "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(laplace3d_neumann_tutorial "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(laplace2d_dirichlet_Ldomain_singularity_tutorial "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(laplace2d_dirichlet_tutorial "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(laplace3d_dirichlet_tutorial "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(ring_of_vorticity "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(yalaplace "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space;toolbox")
# def_demo(laminar_flow_benchmark "applications;dof;fem;linear_algebra;linear_solver;mesh;nonlinear;quadrature;space;timer;toolbox")
# def_demo(instationary_flow_benchmark "applications;dof;fem;linear_algebra;linear_solver;mesh;nonlinear;quadrature;space;toolbox")

# # Stokes demos
# def_demo(stokes2d_convergence_test "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(stokes2d_dirichlet_tutorial "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")
# def_demo(stokes2d_dirichlet_acc "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")

# # Navier-Stokes demos
# def_demo(navier_stokes2d_convergence_test "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space")

# def_demo(2d_lung_laminar "applications;dof;fem;linear_algebra;linear_solver;mesh;nonlinear;quadrature;space;utils;timer;toolbox")

# def_demo(3d_lung_laminar "applications;dof;fem;linear_algebra;linear_solver;mesh;nonlinear;quadrature;space;utils;timer;toolbox;bronchioleTreeModel")

# def_demo(3d_lung_instationary "applications;dof;fem;linear_algebra;linear_solver;mesh;nonlinear;quadrature;space;utils;timer;toolbox")
# def_demo(particle_deposition_test "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space;utils;timer;toolbox")
# def_demo(particle_deposition_test_3d "applications;dof;fem;linear_algebra;linear_solver;mesh;quadrature;space;utils;timer;toolbox")
