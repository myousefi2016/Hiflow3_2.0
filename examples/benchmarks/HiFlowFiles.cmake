set(benchmarks_SOURCES
  assembly_bench
  uniform_refinement_perf_test
)

set(SUBDIRS_BENCHMARKS
  channel_benchmark
  fsi_benchmark
)

if(WITH_HYPRE)
    list(APPEND SUBDIRS_BENCHMARKS channel_benchmark_schur)
endif()

# Recurse through sub-directories and gather information for libhiflow
foreach(dir ${SUBDIRS_BENCHMARKS})
   include(benchmarks/${dir}/HiFlowFiles.cmake)
   make_directory(${PROJECT_BINARY_DIR}/examples/benchmarks/${dir})
   foreach(i ${${dir}_SOURCES})
     add_executable(${i} benchmarks/${dir}/${i}.cc)
     set_target_properties(${i} PROPERTIES RUNTIME_OUTPUT_DIRECTORY benchmarks/${dir})
      target_link_libraries(${i} hiflow)
      if (WITH_HDF5)
        target_link_libraries(${i} hiflow hdf5 dl)
      endif()
      if(WITH_HYPRE)
        target_link_libraries(${i} hiflow HYPRE dl)
      endif()
      if(WITH_HYPRE AND WITH_HDF5)
        target_link_libraries(${i} hiflow HYPRE hdf5 dl)
      endif()
  endforeach()
  foreach(i ${${dir}_ADDITIONAL_FILES})
    file(COPY benchmarks/${dir}/${i}
      DESTINATION ${PROJECT_BINARY_DIR}/examples/benchmarks/${dir})
  endforeach()
  string(TOUPPER "${dir}" DIRECTORY)
endforeach()
