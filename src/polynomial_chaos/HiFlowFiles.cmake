set(polynomial_chaos_SOURCES
  pc_basis.cc
  pc_tensor.cc
  pc_galerkin_vector.cc
  pc_galerkin_matrix.cc
  pc_multilevel.cc
  pc_gmres.cc
)

set(polynomial_chaos_PUBLIC_HEADERS 
  pc_basis.h
  pc_tensor.h
  pc_galerkin_vector.h
  pc_galerkin_matrix.h
  pc_multilevel.h
  pc_residual_assembler.h
  pc_gmres.h
)

if (${WITH_ILUPP})
  list(APPEND polynomial_chaos_SOURCES pc_block_ilupp_preconditioner.cc)
  list(APPEND polynomial_chaos_PUBLIC_HEADERS pc_block_ilupp_preconditioner.h)
endif()

if (${WITH_UMFPACK})
  list(APPEND polynomial_chaos_SOURCES pc_meanbased_preconditioner.cc)
  list(APPEND polynomial_chaos_PUBLIC_HEADERS pc_meanbased_preconditioner.h)
endif()
