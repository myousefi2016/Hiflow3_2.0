include_directories(${PROJECT_BINARY_DIR}/src/include)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR}/contrib/tinyxml)

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

# HDF5
if (WITH_HDF5)
  include_directories(${HDF5_INCLUDE_DIR})
  link_directories(${HDF5_LIBRARY_DIRS})
endif()

if (WITH_CUDA)
  include_directories(${CUDA_TOOLKIT_INCLUDE})
endif()

if (WITH_HYPRE)
  include_directories(${HYPRE_INCLUDE_DIR})
  link_directories(${HYPRE_LIBRARIES})
endif()

if(WITH_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIR})
endif()

if(WITH_PETSC)
  include_directories(${PETSC_INCLUDE_DIR})
endif()

if(WITH_COMPLEX_PETSC)
  include_directories(${PETSC_INCLUDE_DIR})
endif()

if(WITH_SLEPC)
  include_directories(${SLEPC_INCLUDE_DIR})
endif()

if(WITH_P4EST)
  include_directories(${P4EST_INCLUDE_DIR})
  link_directories(${P4EST_SC_LIBRARIES})
  link_directories(${P4EST_Z_LIBRARIES})
  link_directories(${P4EST_LIBRARIES})
endif()

add_definitions(-DMESHES_DATADIR="${PROJECT_BINARY_DIR}/examples/data/")
add_definitions(-DHIFLOW_DIR="${PROJECT_SOURCE_DIR}")

# Helper macro to define a demo, whose source code exists in
# ${DEMO_NAME}.cc. Installs source code in
# share/hiflow3/examples/${SUBDIR} directory.
macro(def_demo DEMO_NAME SUBDIR)
  add_executable(${DEMO_NAME} ${DEMO_NAME}.cc)
  target_link_libraries(${DEMO_NAME} hiflow)
  if (WITH_HDF5)
    target_link_libraries(${DEMO_NAME} hiflow hdf5 dl)
  endif()
  install(FILES ${DEMO_NAME}.cc
    DESTINATION ${HIFLOW_EXAMPLES_DIR}/${SUBDIR})
endmacro(def_demo)

set(SUBDIRS_EXAMPLES
  assembly
  benchmark_local_solver
  benchmarks
  blood_flow
  convection_diffusion
  data
  instationary_flow
  laminar_flow
  laplace
  mesh
  distributed_control_poisson
  poisson
  flow
  flow_dg
  poisson_eigenvalue
  poisson_periodic
  poisson_uncertainty
  poisson_adaptive
  newton
  porous_media
  direct_inverse
  hp_dg_elliptic
  elasticity
  geometric_multigrid
  mixed_precision
  poisson_precond
  met_flow
)

# Recurse through sub-directories and gather information for libhiflow
foreach(dir ${SUBDIRS_EXAMPLES})

  include(${dir}/HiFlowFiles.cmake)

  make_directory(${PROJECT_BINARY_DIR}/examples/${dir})

  foreach(i ${${dir}_SOURCES})
    add_executable(${i} ${dir}/${i}.cc)
    set_target_properties(${i} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${dir})
    target_link_libraries(${i} hiflow)
    if (WITH_HDF5)
      target_link_libraries(${i} hiflow hdf5 dl)
    endif()
  endforeach()

  foreach(i ${${dir}_ADDITIONAL_FILES})
    file(COPY ${dir}/${i}
      DESTINATION ${PROJECT_BINARY_DIR}/examples/${dir})
  endforeach()

  string(TOUPPER "${dir}" DIRECTORY)

endforeach()
