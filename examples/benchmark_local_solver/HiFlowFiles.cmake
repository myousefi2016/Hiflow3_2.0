include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${PROJECT_BINARY_DIR}/src/include)
include_directories(${MPI_INCLUDE_PATH})

set(benchmark_local_solver_SOURCES
  benchmark_local_solver
)
