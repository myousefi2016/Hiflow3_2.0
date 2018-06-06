include_directories(linear_algebra/lmp/opencl)

if (${WITH_OPENCL})

set(opencl_SOURCES
  lmatrix_csr_opencl.cc
  lvector_opencl.cc
  mem_opencl.cc
  opencl_global.cc
  opencl_kernel_mapper.cc
  platform_management_opencl.cc)

set(opencl_PUBLIC_HEADERS 
  mem_opencl.h
  opencl_global.h
  opencl_utils.h
  platform_management_opencl.h
  )

else()
  set(opencl_SOURCES "")
  set(opencl_PUBLIC_HEADERS "")
endif()
