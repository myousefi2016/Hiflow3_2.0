include_directories(linear_algebra/lmp/cuda)

set(cuda_SOURCES
  GPUcublas2_CSR_lMatrix.cc
  GPUcublas2_lVector.cc
  lvector_gpu.cc
  lmatrix_csr_gpu.cc
  lmatrix_coo_gpu.cc
)

if (${WITH_CUDA})
  list(APPEND cuda_SOURCES
    mem_cuda.cu 
    GPUscalar_csr_lmatrix.cu
    GPUscalar_tex_csr_lmatrix.cu 
    GPUblas_lvector.cu 
    platform_management_gpu.cu
    cuda_async_iter.cu
    cuda_gmg.cu
  )
endif()

if (NOT ${WITH_CUDA})
  message("Compiling without CUDA support")
  foreach(i ${cuda_SOURCES})
    set_source_files_properties(linear_algebra/lmp/cuda/${i} PROPERTIES LANGUAGE CXX)
    set_source_files_properties(linear_algebra/lmp/cuda/${i} PROPERTIES COMPILE_FLAGS "-DNOCUDA -I${PROJECT_SOURCE_DIR}/src/linear_algebra/lmp -x c++")
  endforeach()
endif()

# TODO
if (${WITH_CUDA})
set(cuda_PUBLIC_HEADERS
  cuda_async_iter.h
  cuda_gmg.h
)
else()
set(cuda_PUBLIC_HEADERS "")
endif()
