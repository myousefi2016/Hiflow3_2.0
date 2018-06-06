add_subdirectory(AMG)
add_subdirectory(CoupledMatrix_MPI)
add_subdirectory(LinearAlgebra)
add_subdirectory(SeqDenseMatrix)

if (WITH_HYPRE)
    add_subdirectory(Hypre)
endif()

if (WITH_PETSC)
    add_subdirectory(PETSc)
endif()

if (WITH_COMPLEX_PETSC)
    add_subdirectory(PETSc)
endif()

#if (WITH_P4EST)
#    add_subdirectory(P4EST)
#endif()
