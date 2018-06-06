set(linear_algebra_SOURCES
  coupled_matrix.cc
  coupled_vector.cc
  couplings.cc
  hypre_matrix.cc
  hypre_vector.cc
  la_couplings.cc
  platform_vec_mat.cc
  seq_dense_matrix.cc
  seq_dense_matrix_cblas.cc
  seq_dense_matrix_mkl.cc
  pce_matrix.cc
  pce_vector.cc
)

set(linear_algebra_PUBLIC_HEADERS
  coupled_matrix.h
  coupled_vector.h
  couplings.h
  hypre_matrix.h
  hypre_vector.h
  la_couplings.h
  la_descriptor.h
  matrix.h
  platform_vec_mat.h
  seq_dense_matrix.h
  coupled_matrix_factory.h
  coupled_vector_factory.h
  vector.h
  pce_matrix.h
  pce_vector.h
)

if (WITH_PETSC)
    list(APPEND linear_algebra_SOURCES
      petsc_environment.cc
      petsc_matrix.cc
      petsc_vector.cc
    )
    list(APPEND linear_algebra_PUBLIC_HEADERS
      petsc_environment.h
#      petsc_matrix_interface.h NOTE: including these headers results in including petsc.h which results in ambiguous definition of vec (Hiflow vec, PETSc vec). 
#      petsc_vector_interface.h NOTE: therefore: if needed, include headers in correspodning .cc files
      petsc_la_descriptor.h
      petsc_matrix.h
      petsc_vector.h
    )
endif()

# add lmp as include directory
include_directories(linear_algebra/lmp)
include_directories(tools)

# add files from lmp directory
include(linear_algebra/lmp/HiFlowFiles.cmake)

foreach (i ${lmp_SOURCES})
  list(APPEND linear_algebra_SOURCES "lmp/${i}")
endforeach()

foreach (i ${lmp_PUBLIC_HEADERS})
  list(APPEND linear_algebra_PUBLIC_HEADERS "lmp/${i}")
endforeach()
