set(linear_solver_SOURCES
  bicgstab.cc
  cg.cc
  fgmres.cc
  gmres.cc
  hypre_cg.cc
  hypre_bicgstab.cc
  hypre_boomer_amg.cc
  hypre_gmres.cc
  hypre_preconditioner_euclid.cc
  hypre_preconditioner_pilut.cc
  hypre_preconditioner_parasails.cc
  preconditioner_mc.cc
  preconditioner_ai.cc
  preconditioner_bjacobi_standard.cc
  preconditioner_ilupp.cc
  preconditioner_parallel_naive.cc
  preconditioner_vanka.cc
  schur_complement.cc
  mpir.cc
  jacobi.cc
  tri_block.cc
)

set(linear_solver_PUBLIC_HEADERS
  bicgstab.h
  cg.h
  fgmres.h
  gmres.h
  hypre_linear_solver.h
  hypre_preconditioner.h
  hypre_cg.h
  hypre_bicgstab.h
  hypre_boomer_amg.h
  hypre_gmres.h
  hypre_preconditioner_euclid.h
  hypre_preconditioner_pilut.h
  hypre_preconditioner_parasails.h
  linear_solver.h
  preconditioner.h
  preconditioner_bjacobi.h
  preconditioner_bjacobi_standard.h
  preconditioner_ilupp.h
  preconditioner_vanka.h
  preconditioner_mc.h
  preconditioner_ai.h
  linear_solver_creator.h
  linear_solver_factory.h
  preconditioner.h
  preconditioner_parallel_naive.h
  schur_complement.h
  mpir.h
  jacobi.h
  tri_block.h
)

if (WITH_PETSC)
  list(APPEND linear_solver_SOURCES
    petsc_general_ksp.cc
  )
  list(APPEND linear_solver_PUBLIC_HEADERS
    petsc_general_ksp.h
    petsc_linear_solver.h
    petsc_preconditioner.h
  )
endif()

# add files from specific directory
include(linear_solver/gmg/HiFlowFiles.cmake)

foreach (i ${gmg_SOURCES})
  list(APPEND linear_solver_SOURCES "gmg/${i}")
endforeach()

foreach (i ${gmg_PUBLIC_HEADERS})
  list(APPEND linear_solver_PUBLIC_HEADERS "gmg/${i}")
endforeach()

include(linear_solver/gpu-based/HiFlowFiles.cmake)

foreach (i ${gpu-based_SOURCES})
  list(APPEND linear_solver_SOURCES "gpu-based/${i}")
endforeach()

foreach (i ${gpu-based_PUBLIC_HEADERS})
  list(APPEND linear_solver_PUBLIC_HEADERS "gpu-based/${i}")
endforeach()

if (${WITH_MUMPS})
  list(APPEND linear_solver_SOURCES mumps_solver.cc)
  list(APPEND linear_solver_PUBLIC_HEADERS mumps_solver.h mumps_structure.h)
endif()

if (${WITH_ILUPP})
  list(APPEND linear_solver_SOURCES preconditioner_ilupp.cc)
  list(APPEND linear_solver_PUBLIC_HEADERS preconditioner_ilupp.h)
endif()

if (${WITH_UMFPACK})
  list(APPEND linear_solver_SOURCES umfpack_solver.cc)
  list(APPEND linear_solver_PUBLIC_HEADERS umfpack_solver.h)
endif()
